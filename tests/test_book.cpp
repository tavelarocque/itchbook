// ReferenceBook: unit tests of each event's semantics, then a property test
// against a brute-force model that recomputes every level, queue position and
// the checksum from a flat list of orders after every operation.

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/reference_book.hpp"

namespace itchbook::book {
namespace {

std::vector<LevelInfo> L(std::initializer_list<LevelInfo> v) { return v; }

TEST(ReferenceBook, AddBuildsLevelsInPriorityOrder) {
  ReferenceBook b;
  b.add(1, 10, 'B', 100, 50'0000);
  b.add(1, 11, 'B', 200, 51'0000);
  b.add(1, 12, 'B', 300, 50'0000);
  b.add(1, 20, 'S', 400, 52'0000);
  b.add(1, 21, 'S', 500, 53'0000);
  EXPECT_EQ(b.levels(1, Side::buy), L({{51'0000, 200, 1}, {50'0000, 400, 2}}));
  EXPECT_EQ(b.levels(1, Side::sell), L({{52'0000, 400, 1}, {53'0000, 500, 1}}));
  EXPECT_EQ(b.best(1, Side::buy)->price, 51'0000U);
  EXPECT_FALSE(b.best(2, Side::buy).has_value());
  EXPECT_EQ(b.live_orders(), 5U);
  EXPECT_EQ(b.level_count(), 4U);
}

TEST(ReferenceBook, QueuePositionIsTimePriority) {
  ReferenceBook b;
  b.add(1, 1, 'S', 100, 7);
  b.add(1, 2, 'S', 200, 7);
  b.add(1, 3, 'S', 300, 7);
  const auto o = b.order(3);
  ASSERT_TRUE(o);
  EXPECT_EQ(o->rank, 2U);
  EXPECT_EQ(o->shares_ahead, 300U);
  b.cancel(1, 1, 60);  // partial cancel keeps priority
  EXPECT_EQ(b.order(1)->rank, 0U);
  EXPECT_EQ(b.order(3)->shares_ahead, 240U);
}

TEST(ReferenceBook, ExecuteReducesAndRemovesAtZero) {
  ReferenceBook b;
  b.add(1, 1, 'B', 100, 7);
  b.add(1, 2, 'B', 100, 7);
  b.execute(1, 1, 40);
  EXPECT_EQ(b.levels(1, Side::buy), L({{7, 160, 2}}));
  b.execute(1, 1, 60);
  EXPECT_EQ(b.levels(1, Side::buy), L({{7, 100, 1}}));
  EXPECT_FALSE(b.order(1));
  EXPECT_EQ(b.stats().executed_at_head, 2U);
  b.execute(1, 2, 100);
  EXPECT_TRUE(b.levels(1, Side::buy).empty());
  EXPECT_EQ(b.level_count(), 0U);
}

TEST(ReferenceBook, ExecutionBehindHeadIsCounted) {
  ReferenceBook b;
  b.add(1, 1, 'B', 100, 7);
  b.add(1, 2, 'B', 100, 7);
  b.execute(1, 2, 10);
  EXPECT_EQ(b.stats().executed_behind_head, 1U);
  EXPECT_EQ(b.stats().executed_at_head, 0U);
}

TEST(ReferenceBook, ReplaceMovesToTailWithNewRefAndKeepsSide) {
  ReferenceBook b;
  b.add(1, 1, 'S', 100, 7);
  b.add(1, 2, 'S', 100, 7);
  b.replace(1, 1, 99, 150, 7);  // same price, new ref: loses priority
  EXPECT_FALSE(b.order(1));
  const auto o = b.order(99);
  ASSERT_TRUE(o);
  EXPECT_EQ(o->side, Side::sell);
  EXPECT_EQ(o->rank, 1U);
  EXPECT_EQ(o->shares, 150U);  // absolute size, not a delta
  b.replace(1, 99, 100, 10, 8);
  EXPECT_EQ(b.levels(1, Side::sell), L({{7, 100, 1}, {8, 10, 1}}));
}

TEST(ReferenceBook, AnomaliesAreCountedAndLeaveStateUnchanged) {
  ReferenceBook b;
  b.add(1, 1, 'B', 100, 7);
  const auto before = b.checksum();
  b.execute(1, 42, 1);
  b.cancel(1, 42, 1);
  b.remove(1, 42);
  b.replace(1, 42, 43, 1, 1);
  b.add(1, 1, 'B', 5, 5);   // duplicate
  b.add(1, 50, 'Z', 5, 5);  // bad side
  b.add(1, 60, 'S', 5, 9);
  b.replace(1, 1, 60, 5, 5);  // new ref already live
  b.remove(1, 60);
  EXPECT_EQ(b.checksum(), before);
  const auto& s = b.stats();
  EXPECT_EQ(s.orphans(), 4U);
  EXPECT_EQ(s.duplicate_add, 1U);
  EXPECT_EQ(s.bad_side, 1U);
  EXPECT_EQ(s.duplicate_replace, 1U);
}

TEST(ReferenceBook, OverExecuteIsClampedAndCounted) {
  ReferenceBook b;
  b.add(1, 1, 'B', 100, 7);
  b.execute(1, 1, 150);
  EXPECT_EQ(b.stats().over_execute, 1U);
  EXPECT_EQ(b.live_orders(), 0U);
}

TEST(ReferenceBook, LocateMismatchIsCountedButOrderLocateWins) {
  ReferenceBook b;
  b.add(3, 1, 'B', 100, 7);
  b.execute(4, 1, 10);
  EXPECT_EQ(b.stats().locate_mismatch, 1U);
  EXPECT_EQ(b.levels(3, Side::buy), L({{7, 90, 1}}));
}

// ---- brute-force model -----------------------------------------------------

struct ModelOrder {
  Locate locate;
  Side side;
  Price price;
  Shares shares;
  std::uint64_t seq;  // insertion order: lower = earlier in queue
};

struct Model {
  std::map<OrderRef, ModelOrder> orders;
  std::uint64_t seq = 0;
  std::uint64_t orphans = 0;

  void add(Locate l, OrderRef r, char s, Shares n, Price p) {
    if ((s != 'B' && s != 'S') || orders.contains(r)) {
      return;
    }
    orders[r] = {l, s == 'B' ? Side::buy : Side::sell, p, n, seq++};
  }
  void reduce(OrderRef r, Shares n) {
    auto it = orders.find(r);
    if (it == orders.end()) {
      ++orphans;
      return;
    }
    if (n >= it->second.shares) {
      orders.erase(it);
    } else {
      it->second.shares -= n;
    }
  }
  void remove(OrderRef r) {
    if (orders.erase(r) == 0) {
      ++orphans;
    }
  }
  void replace(OrderRef o, OrderRef r, Shares n, Price p) {
    auto it = orders.find(o);
    if (it == orders.end()) {
      ++orphans;
      return;
    }
    if (r != o && orders.contains(r)) {
      return;
    }
    const ModelOrder old = it->second;
    orders.erase(it);
    if (n > 0) {
      orders[r] = {old.locate, old.side, p, n, seq++};
    }
  }

  // Orders of one (locate, side) in book priority: price, then time.
  [[nodiscard]] std::vector<std::pair<OrderRef, ModelOrder>> sorted(Locate l, Side s) const {
    std::vector<std::pair<OrderRef, ModelOrder>> v;
    for (const auto& [r, o] : orders) {
      if (o.locate == l && o.side == s) {
        v.emplace_back(r, o);
      }
    }
    std::sort(v.begin(), v.end(), [&](const auto& a, const auto& b) {
      if (a.second.price != b.second.price) {
        return s == Side::buy ? a.second.price > b.second.price : a.second.price < b.second.price;
      }
      return a.second.seq < b.second.seq;
    });
    return v;
  }
  [[nodiscard]] std::vector<LevelInfo> levels(Locate l, Side s) const {
    std::vector<LevelInfo> out;
    for (const auto& [r, o] : sorted(l, s)) {
      if (out.empty() || out.back().price != o.price) {
        out.push_back({o.price, 0, 0});
      }
      out.back().shares += o.shares;
      ++out.back().orders;
    }
    return out;
  }
  [[nodiscard]] std::uint64_t checksum(Locate max_locate) const {
    std::uint64_t sum = 0;
    for (Locate l = 0; l <= max_locate; ++l) {
      for (const Side s : {Side::buy, Side::sell}) {
        for (const auto& lv : levels(l, s)) {
          sum += mix({l, static_cast<std::uint64_t>(s), lv.price, lv.shares, lv.orders});
        }
        Price cur = 0;
        std::uint64_t rank = 0;
        bool first = true;
        for (const auto& [r, o] : sorted(l, s)) {
          if (first || o.price != cur) {
            rank = 0;
            cur = o.price;
            first = false;
          }
          sum += mix({r, l, static_cast<std::uint64_t>(s), o.price, o.shares, rank++});
        }
      }
    }
    return sum;
  }
};

RC_GTEST_PROP(ReferenceBookProperties, MatchesBruteForceModel, ()) {
  constexpr Locate kMaxLocate = 3;
  const auto n_ops = *rc::gen::inRange(1, 400);
  ReferenceBook b;
  Model m;
  OrderRef next_ref = 1;
  std::vector<OrderRef> issued;

  for (int i = 0; i < n_ops; ++i) {
    const int kind = *rc::gen::inRange(0, 100);
    const auto l = *rc::gen::inRange<Locate>(1, kMaxLocate + 1);
    const auto price = *rc::gen::inRange<Price>(1, 8);  // few prices, so levels hold several orders
    const auto shares = *rc::gen::inRange<Shares>(0, 400);
    // Mostly refer to an order that was issued (live or not); sometimes a fresh unknown one.
    const OrderRef target = (!issued.empty() && *rc::gen::inRange(0, 10) < 9)
                                ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())]
                                : next_ref + 1000;
    const Locate header = m.orders.contains(target) ? m.orders.at(target).locate : l;
    if (kind < 40) {
      const char side = *rc::gen::element('B', 'S', 'B', 'S', 'X');
      b.add(l, next_ref, side, shares == 0 ? 1 : shares, price);
      m.add(l, next_ref, side, shares == 0 ? 1 : shares, price);
      issued.push_back(next_ref++);
    } else if (kind < 55) {
      b.execute(header, target, shares);
      m.reduce(target, shares);
    } else if (kind < 70) {
      b.cancel(header, target, shares);
      m.reduce(target, shares);
    } else if (kind < 85) {
      b.remove(header, target);
      m.remove(target);
    } else {
      b.replace(header, target, next_ref, shares, price);
      m.replace(target, next_ref, shares, price);
      issued.push_back(next_ref++);
    }

    for (Locate ll = 1; ll <= kMaxLocate; ++ll) {
      for (const Side s : {Side::buy, Side::sell}) {
        RC_ASSERT(b.levels(ll, s) == m.levels(ll, s));
      }
    }
    RC_ASSERT(b.live_orders() == m.orders.size());
  }
  RC_ASSERT(b.checksum() == m.checksum(kMaxLocate));
  RC_ASSERT(b.stats().orphans() == m.orphans);
  // Queue position of every live order matches the model's ordering.
  for (const auto& [r, o] : m.orders) {
    const auto info = b.order(r);
    RC_ASSERT(info.has_value());
    std::size_t rank = 0;
    std::uint64_t ahead = 0;
    for (const auto& [r2, o2] : m.sorted(o.locate, o.side)) {
      if (o2.price != o.price) {
        continue;
      }
      if (r2 == r) {
        break;
      }
      ++rank;
      ahead += o2.shares;
    }
    RC_ASSERT(info->rank == rank);
    RC_ASSERT(info->shares_ahead == ahead);
  }
}

}  // namespace
}  // namespace itchbook::book

// ---- FastBook: differential testing against ReferenceBook ----------------------

namespace itchbook::book {
namespace {

void expect_same(const ReferenceBook& r, const FastBook& f, Locate max_locate) {
  for (Locate l = 0; l <= max_locate; ++l) {
    for (const Side s : {Side::buy, Side::sell}) {
      RC_ASSERT(r.levels(l, s) == f.levels(l, s));
      RC_ASSERT(r.best(l, s) == f.best(l, s));
    }
  }
  RC_ASSERT(r.live_orders() == f.live_orders());
  RC_ASSERT(r.level_count() == f.level_count());
}

bool same_stats(const BookStats& a, const BookStats& b) {
  return a.adds == b.adds && a.executes == b.executes && a.cancels == b.cancels && a.deletes == b.deletes &&
         a.replaces == b.replaces && a.orphan_execute == b.orphan_execute && a.orphan_cancel == b.orphan_cancel &&
         a.orphan_delete == b.orphan_delete && a.orphan_replace == b.orphan_replace &&
         a.duplicate_add == b.duplicate_add && a.duplicate_replace == b.duplicate_replace && a.bad_side == b.bad_side &&
         a.over_execute == b.over_execute && a.over_cancel == b.over_cancel && a.locate_mismatch == b.locate_mismatch &&
         a.replace_to_zero == b.replace_to_zero && a.executed_at_head == b.executed_at_head &&
         a.executed_behind_head == b.executed_behind_head && a.cancel_to_zero == b.cancel_to_zero &&
         a.max_live_orders == b.max_live_orders;
}

// Same event stream into both books. `n_prices` large forces the long-book
// binary-search path; a tiny expected size forces the hash index to grow.
void differential(int n_ops, Price n_prices, std::size_t expected) {
  constexpr Locate kMaxLocate = 3;
  ReferenceBook r;
  FastBook f(expected);
  OrderRef next_ref = 1;
  std::vector<OrderRef> issued;
  for (int i = 0; i < n_ops; ++i) {
    const int kind = *rc::gen::inRange(0, 100);
    const auto l = *rc::gen::inRange<Locate>(1, kMaxLocate + 1);
    const auto price = *rc::gen::inRange<Price>(1, n_prices + 1);
    const auto shares = *rc::gen::inRange<Shares>(0, 400);
    const OrderRef target = (!issued.empty() && *rc::gen::inRange(0, 10) < 9)
                                ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())]
                                : next_ref + 1000;
    const auto header = *rc::gen::inRange<Locate>(1, kMaxLocate + 1);
    if (kind < 45) {
      const char side = *rc::gen::element('B', 'S', 'B', 'S', 'X');
      r.add(l, next_ref, side, shares, price);
      f.add(l, next_ref, side, shares, price);
      issued.push_back(next_ref++);
    } else if (kind < 60) {
      RC_ASSERT(r.at_head(target) == f.at_head(target));
      r.execute(header, target, shares);
      f.execute(header, target, shares);
    } else if (kind < 70) {
      r.cancel(header, target, shares);
      f.cancel(header, target, shares);
    } else if (kind < 85) {
      r.remove(header, target);
      f.remove(header, target);
    } else {
      const OrderRef nr = *rc::gen::inRange(0, 20) == 0 && !issued.empty() ? issued.back() : next_ref;
      r.replace(header, target, nr, shares, price);
      f.replace(header, target, nr, shares, price);
      issued.push_back(next_ref++);
    }
    expect_same(r, f, kMaxLocate);
  }
  RC_ASSERT(r.checksum() == f.checksum());
  RC_ASSERT(same_stats(r.stats(), f.stats()));
  for (const OrderRef ref : issued) {
    const auto a = r.order(ref);
    const auto b = f.order(ref);
    RC_ASSERT(a.has_value() == b.has_value());
    if (a) {
      RC_ASSERT(a->rank == b->rank);
      RC_ASSERT(a->shares_ahead == b->shares_ahead);
      RC_ASSERT(a->shares == b->shares);
    }
  }
}

RC_GTEST_PROP(FastBookProperties, MatchesReferenceOnCrowdedLevels, ()) {
  differential(*rc::gen::inRange(1, 500), 6, 1024);
}
RC_GTEST_PROP(FastBookProperties, MatchesReferenceOnDeepBooks, ()) {
  differential(*rc::gen::inRange(1, 800), 200, 1024);
}
RC_GTEST_PROP(FastBookProperties, MatchesReferenceWhileIndexGrows, ()) {
  differential(*rc::gen::inRange(600, 1500), 40, 1);
}

}  // namespace
}  // namespace itchbook::book

// ---- order indexes against std::unordered_map ----------------------------------

namespace itchbook::book {
namespace {

template <class Index>
void index_matches_map() {
  Index idx(16);
  std::unordered_map<OrderRef, std::uint32_t> ref;
  const auto n = *rc::gen::inRange(1, 3000);
  OrderRef next = *rc::gen::inRange<OrderRef>(0, 1U << 20U);
  std::vector<OrderRef> live;
  for (int i = 0; i < n; ++i) {
    const int op = *rc::gen::inRange(0, 10);
    if (op < 5 || live.empty()) {
      // Mostly increasing refs, like NASDAQ; sometimes a jump, sometimes a huge (fallback) ref.
      const int kind = *rc::gen::inRange(0, 20);
      OrderRef r = kind == 0 ? (OrderRef{1} << 40U) + *rc::gen::inRange<OrderRef>(0, 1000)
                             : (kind == 1 ? next + 200'000 : next + *rc::gen::inRange<OrderRef>(1, 50));
      next = std::max(next, r < (OrderRef{1} << 40U) ? r : next);
      if (ref.contains(r)) {
        continue;
      }
      const auto slot = static_cast<std::uint32_t>(i);
      idx.insert(r, slot);
      ref[r] = slot;
      live.push_back(r);
    } else {
      const std::size_t k = *rc::gen::inRange<std::size_t>(0, live.size());
      idx.erase(live[k]);
      ref.erase(live[k]);
      live[k] = live.back();
      live.pop_back();
    }
    const OrderRef probe = *rc::gen::inRange<OrderRef>(0, next + 10);
    RC_ASSERT(idx.find(probe) == (ref.contains(probe) ? ref[probe] : kNoSlot));
  }
  for (const auto& [r, s] : ref) {
    RC_ASSERT(idx.find(r) == s);
  }
}

RC_GTEST_PROP(RefIndexProperties, HashMatchesMap, ()) { index_matches_map<HashRefIndex>(); }
RC_GTEST_PROP(RefIndexProperties, PagedMatchesMap, ()) { index_matches_map<PagedRefIndex>(); }

TEST(PagedRefIndex, ReleasesEmptiedPagesBehindTheFrontier) {
  PagedRefIndex idx;
  for (OrderRef r = 0; r < 4 * PagedRefIndex::kPageSize; ++r) {
    idx.insert(r, static_cast<std::uint32_t>(r));
  }
  const auto full = idx.memory_bytes();
  for (OrderRef r = 0; r < 3 * PagedRefIndex::kPageSize; ++r) {
    idx.erase(r);
  }
  EXPECT_LT(idx.memory_bytes(), full / 2);
  EXPECT_EQ(idx.find(3 * PagedRefIndex::kPageSize + 5), 3 * PagedRefIndex::kPageSize + 5);
  EXPECT_EQ(idx.find(5), kNoSlot);
}

}  // namespace
}  // namespace itchbook::book
