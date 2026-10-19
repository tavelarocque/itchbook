#pragma once
// ReferenceBook: the M2 L3 order book, written for clarity and checkability
// rather than speed. It is the in-repo oracle the optimised M3 engine must
// match event for event.
//
// State: every live order keyed by its order reference number; per symbol
// (stock locate) and side, price levels in priority order, each holding its
// orders in time priority (FIFO). Prices stay in ITCH fixed point.
//
// Semantics, applied identically by the Python reference (tools/itchref/bookref.py):
//   add       side must be B/S, ref must be new; append at the tail of its level.
//   execute   (E and C) shares -= n, remove at 0. C never reprices the order.
//   cancel    (X) shares -= n, keeps priority; remove at 0.
//   remove    (D) remove the order.
//   replace   (U) remove the original, insert the new ref with the original's
//             locate and side at the tail of the new level (priority lost).
// Anomalies are counted, never thrown: a reference to an unknown order (an
// "orphan"), a duplicate add, an execute or cancel larger than the remaining
// shares (clamped), or a header stock locate that differs from the order's.

#include <array>
#include <cstdint>
#include <functional>
#include <list>
#include <map>
#include <optional>
#include <unordered_map>
#include <vector>

#include "itchbook/book/types.hpp"

namespace itchbook::book {

struct BookStats {
  std::uint64_t adds = 0;
  std::uint64_t executes = 0;
  std::uint64_t cancels = 0;
  std::uint64_t deletes = 0;
  std::uint64_t replaces = 0;

  std::uint64_t orphan_execute = 0;
  std::uint64_t orphan_cancel = 0;
  std::uint64_t orphan_delete = 0;
  std::uint64_t orphan_replace = 0;
  std::uint64_t duplicate_add = 0;
  std::uint64_t duplicate_replace = 0;  // U whose new ref is already live
  std::uint64_t bad_side = 0;
  std::uint64_t over_execute = 0;
  std::uint64_t over_cancel = 0;
  std::uint64_t locate_mismatch = 0;
  std::uint64_t replace_to_zero = 0;  // U with 0 shares: original removed, nothing inserted

  std::uint64_t executed_at_head = 0;  // execution hit the order first in its queue
  std::uint64_t executed_behind_head = 0;
  std::uint64_t cancel_to_zero = 0;
  std::uint64_t max_live_orders = 0;

  [[nodiscard]] std::uint64_t orphans() const noexcept {
    return orphan_execute + orphan_cancel + orphan_delete + orphan_replace;
  }
  [[nodiscard]] std::uint64_t anomalies() const noexcept {
    return orphans() + duplicate_add + duplicate_replace + bad_side + over_execute + over_cancel + locate_mismatch;
  }
};

class ReferenceBook {
 public:
  ReferenceBook();

  void add(Locate locate, OrderRef ref, char side, Shares shares, Price price);
  void execute(Locate locate, OrderRef ref, Shares shares);
  void cancel(Locate locate, OrderRef ref, Shares shares);
  void remove(Locate locate, OrderRef ref);
  void replace(Locate locate, OrderRef old_ref, OrderRef new_ref, Shares shares, Price price);

  // ---- queries -------------------------------------------------------------
  [[nodiscard]] std::optional<LevelInfo> best(Locate locate, Side side) const;
  // Up to `depth` levels from the best price outward.
  [[nodiscard]] std::vector<LevelInfo> levels(Locate locate, Side side, std::size_t depth = SIZE_MAX) const;
  [[nodiscard]] std::size_t live_orders() const noexcept { return orders_.size(); }
  [[nodiscard]] std::size_t level_count() const noexcept { return level_count_; }

  struct OrderInfo {
    Locate locate;
    Side side;
    Price price;
    Shares shares;
    std::size_t rank;            // 0 = first in queue at its price
    std::uint64_t shares_ahead;  // displayed shares queued in front of it
  };
  // O(queue length): for tests and spot checks, not the hot path.
  [[nodiscard]] std::optional<OrderInfo> order(OrderRef ref) const;
  // O(1): is this live order first in its queue? (false if unknown)
  [[nodiscard]] bool at_head(OrderRef ref) const;

  // Order-sensitive checksum of the whole market: every level's aggregate and
  // every order's position in its queue. See tools/itchref/bookref.py.
  [[nodiscard]] std::uint64_t checksum() const;

  [[nodiscard]] const BookStats& stats() const noexcept { return stats_; }

  // Visit every non-empty book as (locate, bids best-first, asks best-first).
  template <class F>
  void for_each_symbol(F&& f) const {
    for (std::size_t l = 0; l < books_.size(); ++l) {
      const SymbolBook& b = books_[l];
      if (!b.bids.empty() || !b.asks.empty()) {
        f(static_cast<Locate>(l));
      }
    }
  }

 private:
  struct Level {
    std::uint64_t shares = 0;
    std::list<OrderRef> queue;  // time priority, front = first to fill
  };
  struct Order {
    Locate locate;
    Side side;
    Price price;
    Shares shares;
    std::list<OrderRef>::iterator pos;
  };
  struct SymbolBook {
    std::map<Price, Level, std::greater<>> bids;  // best (highest) first
    std::map<Price, Level, std::less<>> asks;     // best (lowest) first
  };

  template <class F>
  decltype(auto) with_side(SymbolBook& b, Side s, F&& f) {
    return s == Side::buy ? f(b.bids) : f(b.asks);
  }
  template <class F>
  decltype(auto) with_side(const SymbolBook& b, Side s, F&& f) const {
    return s == Side::buy ? f(b.bids) : f(b.asks);
  }

  void insert(OrderRef ref, Locate locate, Side side, Shares shares, Price price);
  void erase(std::unordered_map<OrderRef, Order>::iterator it);
  void reduce(std::unordered_map<OrderRef, Order>::iterator it, Shares n);
  void check_locate(Locate header, const Order& o) noexcept;

  std::vector<SymbolBook> books_;  // indexed by stock locate
  std::unordered_map<OrderRef, Order> orders_;
  std::size_t level_count_ = 0;
  BookStats stats_;
};

}  // namespace itchbook::book
