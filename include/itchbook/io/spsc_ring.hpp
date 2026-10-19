#pragma once
// Single-producer single-consumer ring of byte blocks, and a ByteSource that
// reads another source on a background thread through it.
//
// Used by the benchmark to take decompression off the measured thread: the
// producer inflates gzip into blocks while the consumer parses and updates the
// book, and the consumer's own time (minus time spent waiting for a block) is
// what gets reported.
//
// Lock-free protocol: `head_` is written only by the producer and `tail_` only
// by the consumer. The producer fills slot head % N, then publishes it with a
// release store of head + 1; the consumer's acquire load of head_ therefore
// sees the block's bytes. Symmetrically the consumer releases a slot with a
// release store of tail_, which the producer acquires before reusing it. Each
// index sits on its own cache line so the two threads never write the same line.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <memory>
#include <new>
#include <thread>
#include <vector>

#include "itchbook/io/source.hpp"

namespace itchbook::io {

inline constexpr std::size_t kCacheLine = 64;

class SpscBlockRing {
 public:
  SpscBlockRing(std::size_t blocks, std::size_t block_bytes) : n_(blocks), slots_(blocks) {
    for (auto& s : slots_) {
      s.data.resize(block_bytes);
    }
  }

  // Producer side ------------------------------------------------------------
  // Returns the next free block, spinning (with yields) while the ring is full.
  std::vector<std::uint8_t>& begin_write() {
    const std::uint64_t h = head_.load(std::memory_order_relaxed);
    while (h - tail_.load(std::memory_order_acquire) == n_) {
      std::this_thread::yield();
    }
    return slots_[h % n_].data;
  }
  void commit_write(std::size_t bytes) {
    const std::uint64_t h = head_.load(std::memory_order_relaxed);
    slots_[h % n_].size = bytes;
    head_.store(h + 1, std::memory_order_release);
  }

  // Consumer side ------------------------------------------------------------
  // Returns the oldest filled block, or nullptr if none is ready yet.
  struct Block {
    const std::uint8_t* data;
    std::size_t size;
  };
  [[nodiscard]] bool try_read(Block& out) {
    const std::uint64_t t = tail_.load(std::memory_order_relaxed);
    if (t == head_.load(std::memory_order_acquire)) {
      return false;
    }
    const Slot& s = slots_[t % n_];
    out = {s.data.data(), s.size};
    return true;
  }
  void release_read() { tail_.store(tail_.load(std::memory_order_relaxed) + 1, std::memory_order_release); }

 private:
  struct Slot {
    std::vector<std::uint8_t> data;
    std::size_t size = 0;
  };
  const std::size_t n_;
  std::vector<Slot> slots_;
  alignas(kCacheLine) std::atomic<std::uint64_t> head_{0};
  alignas(kCacheLine) std::atomic<std::uint64_t> tail_{0};
};

// Reads `inner` on a background thread. A zero-size block marks end of input.
class PrefetchSource final : public ByteSource {
 public:
  explicit PrefetchSource(std::unique_ptr<ByteSource> inner, std::size_t blocks = 16,
                          std::size_t block_bytes = std::size_t{1} << 20U)
      : inner_(std::move(inner)), ring_(blocks, block_bytes) {
    producer_ = std::thread([this] { produce(); });
  }
  ~PrefetchSource() override {
    stop_.store(true, std::memory_order_relaxed);
    // Drain so a producer blocked on a full ring can see stop_ and exit.
    SpscBlockRing::Block b{};
    while (producer_.joinable()) {
      if (ring_.try_read(b)) {
        ring_.release_read();
      }
      if (done_.load(std::memory_order_acquire)) {
        producer_.join();
      } else {
        std::this_thread::yield();
      }
    }
  }
  PrefetchSource(const PrefetchSource&) = delete;
  PrefetchSource& operator=(const PrefetchSource&) = delete;
  PrefetchSource(PrefetchSource&&) = delete;
  PrefetchSource& operator=(PrefetchSource&&) = delete;

  std::size_t read(std::uint8_t* dst, std::size_t cap) override {
    std::size_t out = 0;
    while (out < cap && !eof_) {
      if (cur_.data == nullptr) {
        if (!ring_.try_read(cur_)) {
          if (out > 0) {
            return out;  // hand over what we have rather than wait
          }
          const auto t0 = std::chrono::steady_clock::now();
          while (!ring_.try_read(cur_)) {
            std::this_thread::yield();
          }
          wait_ns_ += static_cast<std::uint64_t>(
              std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count());
        }
        pos_ = 0;
        if (cur_.size == 0) {
          eof_ = true;
          ring_.release_read();
          cur_ = {};
          break;
        }
      }
      const std::size_t n = std::min(cap - out, cur_.size - pos_);
      std::memcpy(dst + out, cur_.data + pos_, n);
      out += n;
      pos_ += n;
      if (pos_ == cur_.size) {
        ring_.release_read();
        cur_ = {};
      }
    }
    if (eof_ && error_) {
      std::rethrow_exception(error_);
    }
    return out;
  }
  [[nodiscard]] bool truncated() const noexcept override { return truncated_; }
  // Time the consumer spent blocked waiting for the producer.
  [[nodiscard]] std::uint64_t wait_ns() const noexcept { return wait_ns_; }

 private:
  void produce() {
    try {
      for (;;) {
        if (stop_.load(std::memory_order_relaxed)) {
          break;
        }
        auto& block = ring_.begin_write();
        std::size_t filled = 0;
        while (filled < block.size()) {
          const std::size_t n = inner_->read(block.data() + filled, block.size() - filled);
          if (n == 0) {
            break;
          }
          filled += n;
        }
        if (filled == 0) {
          break;
        }
        ring_.commit_write(filled);
      }
    } catch (...) {
      error_ = std::current_exception();
    }
    truncated_ = inner_->truncated();
    ring_.begin_write();
    ring_.commit_write(0);  // end-of-input marker
    done_.store(true, std::memory_order_release);
  }

  std::unique_ptr<ByteSource> inner_;
  SpscBlockRing ring_;
  std::thread producer_;
  std::atomic<bool> stop_{false};
  std::atomic<bool> done_{false};
  std::exception_ptr error_;  // written by producer before the end marker is published
  bool truncated_ = false;    // likewise
  SpscBlockRing::Block cur_{};
  std::size_t pos_ = 0;
  bool eof_ = false;
  std::uint64_t wait_ns_ = 0;
};

}  // namespace itchbook::io
