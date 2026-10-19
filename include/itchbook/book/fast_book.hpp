#pragma once
// FastBook: the M3 engine. Same interface and the same observable behaviour
// as ReferenceBook (every counter, level, queue position and checksum), built
// for speed:
//
//   * orders live in one pooled array and are addressed by 32-bit index; a
//     free list recycles slots, so the hot path never allocates once the pool
//     has grown to the day's peak;
//   * order reference -> slot index is a paged direct-mapped table
//     (PagedRefIndex; -DITCHBOOK_HASH_INDEX selects the open-addressing hash
//     table instead, for comparison);
//   * price levels live in a pool and every order stores its level's index,
//     so executions, cancels and deletes never search for their price; each
//     side of a symbol's book is a vector of (price, level) sorted with the
//     BEST price at the BACK, searched only when a level is created or
//     removed, where inserting and erasing moves the fewest elements;
//   * each level's time-priority queue is an intrusive doubly linked list
//     threaded through the order slots (head, tail, prev, next are indices).
//
// ReferenceBook is the oracle: tests and the full-day gate require identical
// output from both.

#include <algorithm>
#include <cstdint>
#include <optional>
#include <vector>

#include "itchbook/book/ref_index.hpp"
#include "itchbook/book/reference_book.hpp"
#include "itchbook/book/types.hpp"

namespace itchbook::book {

class FastBook {
 public:
#if defined(ITCHBOOK_HASH_INDEX)
  using Index = HashRefIndex;
#else
  using Index = PagedRefIndex;
#endif
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
  [[nodiscard]] const Index& index() const noexcept { return index_; }

  // Software prefetch for an ITCH message that will be applied soon (see
  // itch::Prefetching). Stage 0: the index entry of every order reference the
  // message names, and for adds the symbol's book header. Stage 1: the order
  // slot the index entry points to, and for adds the best price level of the
  // symbol's side. Reads only bytes within `len`.
  void prefetch(const std::uint8_t* p, std::size_t len, int stage) const noexcept;

 private:
  static constexpr std::uint32_t kNil = UINT32_MAX;

  struct Slot {
    OrderRef ref;
    Price price;
    Shares shares;
    std::uint32_t prev;
    std::uint32_t next;   // also the free-list link when the slot is free
    std::uint32_t level;  // index into levels_: no price search on execute/cancel/delete
    Locate locate;
    Side side;
  };
  // Price levels live in a pool and never move, so orders can point at them.
  struct Level {
    std::uint64_t shares;
    Price price;
    std::uint32_t count;
    std::uint32_t head;  // also the free-list link when the level is free
    std::uint32_t tail;
  };
  // Each side of a symbol's book: (price, level) sorted with the BEST price
  // at the back. Only creating or removing a price level searches it.
  struct PriceRef {
    Price price;
    std::uint32_t level;
  };
  using Side_ = std::vector<PriceRef>;
  struct SymbolBook {
    Side_ sides[2];
  };

  // ---- order index ----------------------------------------------------------
  [[nodiscard]] std::uint32_t find(OrderRef ref) const noexcept { return index_.find(ref); }

  // ---- levels ------------------------------------------------------------
  // "Better" means closer to the back of the vector.
  [[nodiscard]] static bool better(Side s, Price a, Price b) noexcept { return s == Side::buy ? a > b : a < b; }
  [[nodiscard]] std::size_t position(const Side_& v, Side s, Price p) const noexcept;
  std::uint32_t level_for_insert(Locate l, Side s, Price p);
  void remove_level(Locate l, Side s, Price p, std::uint32_t level);

  std::uint32_t alloc_slot();
  void insert(OrderRef ref, Locate locate, Side side, Shares shares, Price price);
  void erase(std::uint32_t slot);
  void reduce(std::uint32_t slot, Shares n);

  std::vector<SymbolBook> books_;
  std::vector<Slot> slots_;
  std::uint32_t free_ = kNil;
  std::vector<Level> levels_;
  std::uint32_t free_level_ = kNil;
  Index index_;
  std::size_t live_ = 0;
  std::size_t level_count_ = 0;
  BookStats stats_;
};

}  // namespace itchbook::book
