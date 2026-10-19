#pragma once
// Order reference -> slot index maps for FastBook.
//
// HashRefIndex: open addressing with linear probing and backward-shift
// deletion. Memory is proportional to live orders; every lookup is a random
// access into a table of (2 x live) buckets.
//
// PagedRefIndex: a direct-mapped table split into 64K-entry pages, allocated
// on first use. NASDAQ order references are dense and increase through the
// day (01302020 uses 494M of them), so consecutive adds land in the same
// page and cache line, and an execute or cancel usually refers to a recently
// added order whose page is still in cache. Pages behind the newest one are
// released when their last order dies, so memory follows the live span rather
// than the largest reference. References beyond 2^36 (never seen in real
// files; possible in fuzzed input) go to a HashRefIndex fallback, so a hostile
// reference cannot make the page table explode.

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "itchbook/book/types.hpp"

namespace itchbook::book {

inline constexpr std::uint32_t kNoSlot = UINT32_MAX;

class HashRefIndex {
 public:
  explicit HashRefIndex(std::size_t expected = std::size_t{1} << 20U) {
    const std::size_t size = std::bit_ceil(std::max<std::size_t>(expected * 2, 1024));
    buckets_.assign(size, Bucket{0, kNoSlot});
    mask_ = size - 1;
    shift_ = static_cast<unsigned>(64 - std::countr_zero(size));
  }

  [[nodiscard]] std::uint32_t find(OrderRef ref) const noexcept {
    for (std::size_t i = home(ref);; i = (i + 1) & mask_) {
      const Bucket& b = buckets_[i];
      if (b.slot == kNoSlot) {
        return kNoSlot;
      }
      if (b.ref == ref) {
        return b.slot;
      }
    }
  }
  void insert(OrderRef ref, std::uint32_t slot) {
    if ((size_ + 1) * 2 > buckets_.size()) {
      grow();
    }
    std::size_t i = home(ref);
    while (buckets_[i].slot != kNoSlot) {
      i = (i + 1) & mask_;
    }
    buckets_[i] = {ref, slot};
    ++size_;
  }
  void erase(OrderRef ref) noexcept {
    std::size_t i = home(ref);
    while (buckets_[i].slot == kNoSlot || buckets_[i].ref != ref) {
      i = (i + 1) & mask_;
    }
    // Backward shift: pull later entries of the probe run into the hole when
    // that brings them no further from their home bucket.
    for (std::size_t j = (i + 1) & mask_; buckets_[j].slot != kNoSlot; j = (j + 1) & mask_) {
      if (((j - home(buckets_[j].ref)) & mask_) >= ((j - i) & mask_)) {
        buckets_[i] = buckets_[j];
        i = j;
      }
    }
    buckets_[i].slot = kNoSlot;
    --size_;
  }
  void prefetch(OrderRef ref) const noexcept { __builtin_prefetch(&buckets_[home(ref)]); }
  [[nodiscard]] std::size_t memory_bytes() const noexcept { return buckets_.size() * sizeof(Bucket); }
  [[nodiscard]] static constexpr const char* name() noexcept { return "hash"; }

 private:
  struct Bucket {
    OrderRef ref;
    std::uint32_t slot;
  };
  [[nodiscard]] std::size_t home(OrderRef ref) const noexcept {
    return static_cast<std::size_t>((ref * 0x9E3779B97F4A7C15ULL) >> shift_);
  }
  void grow() {
    const std::vector<Bucket> old = std::move(buckets_);
    buckets_.assign(old.size() * 2, Bucket{0, kNoSlot});
    mask_ = buckets_.size() - 1;
    --shift_;
    for (const Bucket& b : old) {
      if (b.slot != kNoSlot) {
        std::size_t i = home(b.ref);
        while (buckets_[i].slot != kNoSlot) {
          i = (i + 1) & mask_;
        }
        buckets_[i] = b;
      }
    }
  }
  std::vector<Bucket> buckets_;
  std::size_t mask_ = 0;
  unsigned shift_ = 0;
  std::size_t size_ = 0;
};

class PagedRefIndex {
 public:
  static constexpr unsigned kPageBits = 16;
  static constexpr std::size_t kPageSize = std::size_t{1} << kPageBits;
  static constexpr std::size_t kMaxPages = std::size_t{1} << 20U;  // direct for refs < 2^36

  explicit PagedRefIndex(std::size_t /*expected*/ = 0) : fallback_(1024) {}

  [[nodiscard]] std::uint32_t find(OrderRef ref) const noexcept {
    const std::uint64_t p = ref >> kPageBits;
    if (p < pages_.size()) {
      const std::uint32_t* pg = pages_[p].get();
      return pg != nullptr ? pg[ref & (kPageSize - 1)] : kNoSlot;
    }
    return p >= kMaxPages ? fallback_.find(ref) : kNoSlot;
  }
  void insert(OrderRef ref, std::uint32_t slot) {
    const std::uint64_t p = ref >> kPageBits;
    if (p >= kMaxPages) {
      fallback_.insert(ref, slot);
      return;
    }
    if (p >= pages_.size()) {
      pages_.resize(p + 1);
      live_.resize(p + 1, 0);
    }
    if (!pages_[p]) {
      pages_[p] = std::make_unique_for_overwrite<std::uint32_t[]>(kPageSize);
      std::fill_n(pages_[p].get(), kPageSize, kNoSlot);
      peak_pages_ = std::max(peak_pages_, ++pages_alive_);
    }
    newest_ = std::max<std::size_t>(newest_, p);
    pages_[p][ref & (kPageSize - 1)] = slot;
    ++live_[p];
  }
  void erase(OrderRef ref) noexcept {
    const std::uint64_t p = ref >> kPageBits;
    if (p >= kMaxPages) {
      fallback_.erase(ref);
      return;
    }
    pages_[p][ref & (kPageSize - 1)] = kNoSlot;
    // Release emptied pages behind the newest one. The newest page is kept so
    // a page at the frontier is not freed and reallocated over and over.
    if (--live_[p] == 0 && p < newest_) {
      pages_[p].reset();
      --pages_alive_;
    }
  }
  void prefetch(OrderRef ref) const noexcept {
    const std::uint64_t p = ref >> kPageBits;
    if (p < pages_.size() && pages_[p]) {
      __builtin_prefetch(&pages_[p][ref & (kPageSize - 1)]);
    }
  }
  [[nodiscard]] std::size_t memory_bytes() const noexcept {
    return pages_alive_ * kPageSize * sizeof(std::uint32_t) + pages_.capacity() * sizeof(void*) +
           fallback_.memory_bytes();
  }
  [[nodiscard]] std::size_t peak_memory_bytes() const noexcept {
    return peak_pages_ * kPageSize * sizeof(std::uint32_t) + pages_.capacity() * sizeof(void*);
  }
  [[nodiscard]] static constexpr const char* name() noexcept { return "paged"; }

 private:
  std::vector<std::unique_ptr<std::uint32_t[]>> pages_;
  std::vector<std::uint32_t> live_;
  std::size_t pages_alive_ = 0;
  std::size_t peak_pages_ = 0;
  std::size_t newest_ = 0;
  HashRefIndex fallback_;
};

}  // namespace itchbook::book
