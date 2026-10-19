#pragma once
// FastBook: the M3 engine. Same interface and the same observable behaviour
// as ReferenceBook (every counter, level, queue position and checksum), built
// for speed:
//
//   * orders live in one pooled array and are addressed by 32-bit index; a
//     free list recycles slots, so the hot path never allocates once the pool
//     has grown to the day's peak;
//   * order reference -> slot index is an open-addressing hash table with
//     linear probing and backward-shift deletion (no tombstones);
//   * each side of each symbol's book is a vector of price levels sorted with
//     the BEST price at the BACK, so the busiest levels sit where inserting
//     and erasing moves the fewest elements;
//   * each level's time-priority queue is an intrusive doubly linked list
//     threaded through the order slots (head, tail, prev, next are indices).
//
// ReferenceBook is the oracle: tests and the full-day gate require identical
// output from both.

#include <algorithm>
#include <cstdint>
#include <optional>
#include <vector>

#include "itchbook/book/reference_book.hpp"
#include "itchbook/book/types.hpp"

namespace itchbook::book {

class FastBook {
 public:
  using OrderInfo = ReferenceBook::OrderInfo;

  explicit FastBook(std::size_t expected_orders = std::size_t{1} << 21U);

  void add(Locate locate, OrderRef ref, char side, Shares shares, Price price);
  void execute(Locate locate, OrderRef ref, Shares shares);
  void cancel(Locate locate, OrderRef ref, Shares shares);
  void remove(Locate locate, OrderRef ref);
  void replace(Locate locate, OrderRef old_ref, OrderRef new_ref, Shares shares, Price price);

  [[nodiscard]] std::optional<LevelInfo> best(Locate locate, Side side) const;
  [[nodiscard]] std::vector<LevelInfo> levels(Locate locate, Side side, std::size_t depth = SIZE_MAX) const;
  [[nodiscard]] std::size_t live_orders() const noexcept { return live_; }
  [[nodiscard]] std::size_t level_count() const noexcept { return level_count_; }
  [[nodiscard]] std::optional<OrderInfo> order(OrderRef ref) const;
  [[nodiscard]] bool at_head(OrderRef ref) const;
  [[nodiscard]] std::uint64_t checksum() const;
  [[nodiscard]] const BookStats& stats() const noexcept { return stats_; }

  // Software prefetch for an ITCH message that will be applied soon (see
  // itch::Prefetching). Stage 0: the hash bucket of every order reference the
  // message names. Stage 1: the order slot that bucket points to, and for adds
  // the best level of the symbol's side. Reads only bytes within `len`.
  void prefetch(const std::uint8_t* p, std::size_t len, int stage) const noexcept;

 private:
  static constexpr std::uint32_t kNil = UINT32_MAX;

  struct Slot {
    OrderRef ref;
    Price price;
    Shares shares;
    std::uint32_t prev;
    std::uint32_t next;  // also the free-list link when the slot is free
    Locate locate;
    Side side;
  };
  struct Level {
    Price price;
    std::uint32_t count;
    std::uint64_t shares;
    std::uint32_t head;
    std::uint32_t tail;
  };
  using Side_ = std::vector<Level>;  // sorted, best price at back()
  struct SymbolBook {
    Side_ sides[2];
  };

  // ---- order index: open addressing ---------------------------------------
  struct Bucket {
    OrderRef ref;
    std::uint32_t slot;  // kNil = empty
  };
  [[nodiscard]] std::size_t home(OrderRef ref) const noexcept {
    return static_cast<std::size_t>((ref * 0x9E3779B97F4A7C15ULL) >> shift_);
  }
  [[nodiscard]] std::uint32_t find(OrderRef ref) const noexcept;
  void index_insert(OrderRef ref, std::uint32_t slot);
  void index_erase(OrderRef ref);
  void grow_index();

  // ---- levels ------------------------------------------------------------
  // "Better" means closer to the back of the vector.
  [[nodiscard]] static bool better(Side s, Price a, Price b) noexcept { return s == Side::buy ? a > b : a < b; }
  [[nodiscard]] Level* find_level(Locate l, Side s, Price p) noexcept;
  [[nodiscard]] const Level* find_level(Locate l, Side s, Price p) const noexcept;
  Level& level_for_insert(Locate l, Side s, Price p);

  std::uint32_t alloc_slot();
  void insert(OrderRef ref, Locate locate, Side side, Shares shares, Price price);
  void erase(std::uint32_t slot);
  void reduce(std::uint32_t slot, Shares n);

  std::vector<SymbolBook> books_;
  std::vector<Slot> slots_;
  std::uint32_t free_ = kNil;
  std::vector<Bucket> index_;
  std::size_t mask_ = 0;
  unsigned shift_ = 0;
  std::size_t live_ = 0;
  std::size_t level_count_ = 0;
  BookStats stats_;
};

}  // namespace itchbook::book
