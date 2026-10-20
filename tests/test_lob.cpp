#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include "itchbook/book/reference_book.hpp"
#include "itchbook/lob/recorder.hpp"

namespace itchbook::lob {
namespace {

constexpr std::uint32_t kInf = std::numeric_limits<std::uint32_t>::max();

Top side(std::initializer_list<std::pair<std::uint32_t, std::int64_t>> lv, bool ask) {
  Top t;
  if (ask) {
    t.px.fill(kInf);
  }
  std::size_t m = 0;
  for (const auto& [p, q] : lv) {
    t.px[m] = p;
    t.sz[m++] = q;
  }
  return t;
}

// Cont, Kukanov & Stoikov (2014), eq. (1), at the best level.
TEST(Ofi, BestLevelCases) {
  const Top ask = side({{1010, 50}}, true);
  std::array<std::int64_t, kDepth> e{};
  ofi_update(side({{1000, 30}}, false), side({{1000, 45}}, false), ask, ask, e);  // bid size up at same price
  EXPECT_EQ(e[0], 15);
  e.fill(0);
  ofi_update(side({{1000, 30}}, false), side({{1000, 10}}, false), ask, ask, e);  // bid size down
  EXPECT_EQ(e[0], -20);
  e.fill(0);
  ofi_update(side({{1000, 30}}, false), side({{1005, 7}}, false), ask, ask, e);  // bid price improves: +new
  EXPECT_EQ(e[0], 7);
  e.fill(0);
  ofi_update(side({{1000, 30}}, false), side({{990, 80}}, false), ask, ask, e);  // bid level gone: -old
  EXPECT_EQ(e[0], -30);
  e.fill(0);
  const Top bid = side({{1000, 30}}, false);
  ofi_update(bid, bid, side({{1010, 50}}, true), side({{1010, 20}}, true), e);  // ask size down: buying pressure
  EXPECT_EQ(e[0], 30);
  e.fill(0);
  ofi_update(bid, bid, side({{1010, 50}}, true), side({{1009, 5}}, true), e);  // ask price improves: -new
  EXPECT_EQ(e[0], -5);
  e.fill(0);
  ofi_update(bid, bid, side({}, true), side({{1010, 40}}, true), e);  // ask level appears from nothing
  EXPECT_EQ(e[0], -40);
}

RecorderOptions opts(std::vector<std::string> syms, std::uint64_t grid = 100, std::uint64_t start = 1000,
                     std::uint64_t end = 2000) {
  RecorderOptions o;
  o.symbols = std::move(syms);
  o.grid_ns = grid;
  o.start_ns = start;
  o.end_ns = end;
  return o;
}

TEST(LobRecorder, GridConventionAndFlows) {
  LobRecorder r(opts({"AAA"}, 100, 1000, 1300));  // three rows: (1000,1100], (1100,1200], (1200,1300]
  ASSERT_EQ(r.rows(), 3u);
  r.directory(1, "AAA");
  r.trading_state(1, 'T');
  auto at = [&](std::uint64_t ts, auto f) {
    r.clock(ts);
    f(ts);
  };
  at(500, [&](auto ts) { r.add(ts, 1, 1, 'B', 100, 5000); });  // before the grid: book only
  at(900, [&](auto ts) { r.add(ts, 1, 2, 'S', 70, 5100); });
  at(1100, [&](auto ts) { r.add(ts, 1, 3, 'B', 40, 5000); });  // exactly t_1: belongs to row 1
  at(1150, [&](auto ts) { r.execute(ts, 1, 2, 30); });         // resting sell hit: buyer-initiated
  at(1160, [&](auto ts) { r.hidden_trade(ts, 1, 'S', 9); });
  at(1250, [&](auto ts) { r.execute(ts, 1, 1, 100); });  // resting buy hit: seller-initiated
  r.finish();
  const auto& rows = r.rows_of(0);
  ASSERT_EQ(rows.size(), 3u);
  EXPECT_EQ(rows[0][kBidPx], 5000);
  EXPECT_EQ(rows[0][kBidSz], 140);  // includes the add at exactly 1100
  EXPECT_EQ(rows[0][kOfi], 40);
  EXPECT_EQ(rows[0][kEvents], 1);
  EXPECT_EQ(rows[1][kAskSz], 40);
  EXPECT_EQ(rows[1][kBuyVol], 30);
  EXPECT_EQ(rows[1][kOfi], 30);  // ask size fell by 30
  EXPECT_EQ(rows[1][kHiddenS], 9);
  EXPECT_EQ(rows[2][kSellVol], 100);
  EXPECT_EQ(rows[2][kBidSz], 40);
  EXPECT_EQ(rows[2][kOfi], -100);
  EXPECT_EQ(rows[2][kTrading], 1);
}

// Brute force: replay the same events into ReferenceBook and recompute every
// row from scratch: top levels at each grid time, OFI from the book's levels
// before and after each event, and signed volume from the resting order's side.
struct Brute {
  book::ReferenceBook b;
  std::map<std::uint64_t, std::array<std::int64_t, kDepth>> ofi;  // row index -> accumulated OFI
  std::map<std::uint64_t, std::pair<std::int64_t, std::int64_t>> vol;

  static Top top(const book::ReferenceBook& b, Locate l, Side s) {
    Top t;
    if (s == Side::sell) {
      t.px.fill(kInf);
    }
    const auto lv = b.levels(l, s, kDepth);
    for (std::size_t m = 0; m < lv.size(); ++m) {
      t.px[m] = lv[m].price;
      t.sz[m] = static_cast<std::int64_t>(lv[m].shares);
    }
    return t;
  }
  template <class F>
  void event(std::uint64_t row, Locate l, F&& apply) {
    const Top b0 = top(b, l, Side::buy);
    const Top a0 = top(b, l, Side::sell);
    apply();
    const Top b1 = top(b, l, Side::buy);
    const Top a1 = top(b, l, Side::sell);
    auto& acc = ofi[row];
    for (std::size_t m = 0; m < kDepth; ++m) {
      // Written out case by case, independently of ofi_update.
      std::int64_t bid_part = 0;
      if (b1.px[m] > b0.px[m]) {
        bid_part = b1.sz[m];
      } else if (b1.px[m] == b0.px[m]) {
        bid_part = b1.sz[m] - b0.sz[m];
      } else {
        bid_part = -b0.sz[m];
      }
      std::int64_t ask_part = 0;
      if (a1.px[m] < a0.px[m]) {
        ask_part = a1.sz[m];
      } else if (a1.px[m] == a0.px[m]) {
        ask_part = a1.sz[m] - a0.sz[m];
      } else {
        ask_part = -a0.sz[m];
      }
      acc[m] += bid_part - ask_part;
    }
  }
};

RC_GTEST_PROP(LobRecorderProperties, MatchesBruteForceReplay, ()) {
  constexpr std::uint64_t kGrid = 50;
  constexpr std::uint64_t kStart = 1000;
  constexpr std::uint64_t kEnd = 3000;
  LobRecorder r(opts({"AAA", "BBB"}, kGrid, kStart, kEnd));
  r.directory(1, "AAA");
  r.directory(2, "BBB");
  r.directory(3, "ZZZ");  // untracked
  Brute br;
  std::vector<OrderRef> issued;
  OrderRef next = 1;
  std::uint64_t ts = 900;
  const auto row_of = [&](std::uint64_t t) -> std::uint64_t {
    return t <= kStart ? 0 : (t - kStart + kGrid - 1) / kGrid;  // 1-based row containing t; 0 = before the grid
  };
  std::vector<std::array<Top, 4>> tops_at_row(r.rows() + 1);  // [row] -> {bidA, askA, bidB, askB}
  std::uint64_t filled_row = 0;
  auto fill_rows = [&](std::uint64_t until_ts) {
    // Rows whose time t_k < until_ts are final once a message at until_ts arrives.
    while (filled_row < r.rows() && kStart + (filled_row + 1) * kGrid < until_ts) {
      ++filled_row;
      tops_at_row[filled_row] = {Brute::top(br.b, 1, Side::buy), Brute::top(br.b, 1, Side::sell),
                                 Brute::top(br.b, 2, Side::buy), Brute::top(br.b, 2, Side::sell)};
    }
  };
  const int n = *rc::gen::inRange(1, 500);
  for (int i = 0; i < n; ++i) {
    ts += *rc::gen::inRange<std::uint64_t>(0, 12);
    fill_rows(ts);
    r.clock(ts);
    const auto l = *rc::gen::inRange<Locate>(1, 4);
    const auto price = *rc::gen::inRange<Price>(100, 106);
    const auto shares = *rc::gen::inRange<Shares>(1, 200);
    const int kind = *rc::gen::inRange(0, 100);
    const OrderRef target = !issued.empty() ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())] : next + 999;
    const auto info = br.b.order(target);
    const Locate tl = info ? info->locate : l;
    const std::uint64_t row = row_of(ts);
    const bool in_grid = ts > kStart && ts <= kEnd;
    const bool tracked = tl == 1 || tl == 2;
    auto ev = [&](Locate loc, auto&& apply) {
      if (in_grid && (loc == 1 || loc == 2)) {
        br.event(row * 4 + loc, loc, apply);
      } else {
        apply();
      }
    };
    if (kind < 45) {
      const char sd = *rc::gen::element('B', 'S');
      r.add(ts, l, next, sd, shares, price);
      ev(l, [&] { br.b.add(l, next, sd, shares, price); });
      issued.push_back(next++);
    } else if (kind < 65) {
      r.execute(ts, l, target, shares);
      if (info && in_grid && tracked) {
        auto& v = br.vol[row * 4 + tl];
        const std::int64_t x = std::min(shares, info->shares);
        (info->side == Side::buy ? v.second : v.first) += x;
      }
      ev(tl, [&] { br.b.execute(l, target, shares); });
    } else if (kind < 75) {
      r.cancel(ts, l, target, shares);
      ev(tl, [&] { br.b.cancel(l, target, shares); });
    } else if (kind < 88) {
      r.remove(ts, l, target);
      ev(tl, [&] { br.b.remove(l, target); });
    } else {
      const auto ns = *rc::gen::inRange<Shares>(0, 200);
      r.replace(ts, l, target, next, ns, price);
      ev(tl, [&] { br.b.replace(l, target, next, ns, price); });
      issued.push_back(next++);
    }
  }
  fill_rows(kEnd + 1);
  r.finish();
  for (std::size_t si = 0; si < 2; ++si) {
    const Locate loc = static_cast<Locate>(si + 1);
    const auto& rows = r.rows_of(si);
    RC_ASSERT(rows.size() == r.rows());
    for (std::uint64_t k = 1; k <= r.rows(); ++k) {
      const Row& row = rows[k - 1];
      const Top& bid = tops_at_row[k][si * 2];
      const Top& ask = tops_at_row[k][si * 2 + 1];
      for (std::size_t m = 0; m < kDepth; ++m) {
        RC_ASSERT(row[kBidPx + m] == static_cast<std::int32_t>(bid.px[m]));
        RC_ASSERT(row[kBidSz + m] == bid.sz[m]);
        RC_ASSERT(row[kAskPx + m] == (ask.px[m] == kInf ? 0 : static_cast<std::int32_t>(ask.px[m])));
        RC_ASSERT(row[kAskSz + m] == ask.sz[m]);
      }
      const auto it = br.ofi.find(k * 4 + loc);
      for (std::size_t m = 0; m < kDepth; ++m) {
        RC_ASSERT(row[kOfi + m] == (it == br.ofi.end() ? 0 : it->second[m]));
      }
      const auto vt = br.vol.find(k * 4 + loc);
      RC_ASSERT(row[kBuyVol] == (vt == br.vol.end() ? 0 : vt->second.first));
      RC_ASSERT(row[kSellVol] == (vt == br.vol.end() ? 0 : vt->second.second));
    }
  }
}

}  // namespace
}  // namespace itchbook::lob
