#pragma once
// LobRecorder: fixed-grid limit-order-book features for a set of symbols, the
// input to the M5 forecasting study.
//
// The grid covers regular hours, 09:30:00 to 16:00:00, in steps of `grid_ns`.
// Row k describes the interval (t_{k-1}, t_k], t_k = start + k * grid_ns:
//
//   * the top kDepth levels of each side as they stand after every message
//     with timestamp <= t_k (the same convention as itch-book's snapshots);
//   * multi-level order-flow imbalance accumulated over every book event of
//     the symbol in the interval, not just between snapshots: for level m,
//
//       e_m = q^b_m 1{P^b_m >= P^b_m'} - q^b_m' 1{P^b_m <= P^b_m'}
//           - q^a_m 1{P^a_m <= P^a_m'} + q^a_m' 1{P^a_m >= P^a_m'}
//
//     where ' marks the level before the event (Cont, Kukanov & Stoikov 2014
//     for m = 1; Xu, Gould & Howison 2018 for deeper levels). A missing bid
//     level has price 0 and a missing ask level an infinite price, both with
//     size 0, so a level appearing or vanishing counts as an arrival or a
//     departure;
//   * traded volume by aggressor side: an execution against a resting buy
//     order is seller-initiated. Hidden executions (P) are kept separately by
//     the side the feed reports;
//   * the number of book events and whether the symbol was trading (state T)
//     at t_k.

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/types.hpp"

namespace itchbook::lob {

using book::Locate;
using book::OrderRef;
using book::Price;
using book::Shares;
using book::Side;

inline constexpr std::size_t kDepth = 10;
// Column layout of a row (all int32).
inline constexpr std::size_t kBidPx = 0;
inline constexpr std::size_t kBidSz = kBidPx + kDepth;
inline constexpr std::size_t kAskPx = kBidSz + kDepth;
inline constexpr std::size_t kAskSz = kAskPx + kDepth;
inline constexpr std::size_t kOfi = kAskSz + kDepth;
inline constexpr std::size_t kBuyVol = kOfi + kDepth;  // buyer-initiated displayed executions
inline constexpr std::size_t kSellVol = kBuyVol + 1;   // seller-initiated displayed executions
inline constexpr std::size_t kHiddenB = kSellVol + 1;  // hidden executions reported with side B
inline constexpr std::size_t kHiddenS = kHiddenB + 1;  // hidden executions reported with side S
inline constexpr std::size_t kEvents = kHiddenS + 1;   // book events in the interval
inline constexpr std::size_t kTrading = kEvents + 1;   // 1 if the symbol was in state T at t_k
inline constexpr std::size_t kCols = kTrading + 1;
using Row = std::array<std::int32_t, kCols>;

inline constexpr std::uint64_t kNsPerSec = 1'000'000'000ULL;
inline constexpr std::uint64_t kOpenNs = 34'200 * kNsPerSec;   // 09:30:00
inline constexpr std::uint64_t kCloseNs = 57'600 * kNsPerSec;  // 16:00:00

struct RecorderOptions {
  std::vector<std::string> symbols;
  std::uint64_t grid_ns = 250'000'000;  // 250 ms
  std::uint64_t start_ns = kOpenNs;
  std::uint64_t end_ns = kCloseNs;
};

// Top of one side: kDepth (price, size) pairs, best first.
struct Top {
  std::array<std::uint32_t, kDepth> px{};
  std::array<std::int64_t, kDepth> sz{};
};

// One event's contribution to multi-level OFI, from the tops before and after.
// Missing bid levels have price 0; missing ask levels must be passed with
// price UINT32_MAX (see effective_ask).
inline void ofi_update(const Top& bid0, const Top& bid1, const Top& ask0, const Top& ask1,
                       std::array<std::int64_t, kDepth>& acc) noexcept {
  for (std::size_t m = 0; m < kDepth; ++m) {
    std::int64_t e = 0;
    if (bid1.px[m] >= bid0.px[m]) {
      e += bid1.sz[m];
    }
    if (bid1.px[m] <= bid0.px[m]) {
      e -= bid0.sz[m];
    }
    if (ask1.px[m] <= ask0.px[m]) {
      e -= ask1.sz[m];
    }
    if (ask1.px[m] >= ask0.px[m]) {
      e += ask0.sz[m];
    }
    acc[m] += e;
  }
}

class LobRecorder {
 public:
  explicit LobRecorder(RecorderOptions o) : o_(std::move(o)) {
    n_rows_ = (o_.end_ns - o_.start_ns) / o_.grid_ns;
    slot_of_.fill(-1);
    syms_.resize(o_.symbols.size());
    for (std::size_t i = 0; i < o_.symbols.size(); ++i) {
      want_.emplace(o_.symbols[i], i);
      syms_[i].rows.reserve(n_rows_);
    }
  }

  // Message clock: call before applying every message, in stream order.
  void clock(std::uint64_t ts) {
    while (next_ <= n_rows_ && ts > grid_time(next_)) {
      emit();
    }
  }
  // After the last message: complete the grid.
  void finish() {
    while (next_ <= n_rows_) {
      emit();
    }
  }

  void directory(Locate l, std::string_view symbol) {
    const auto it = want_.find(std::string(symbol));
    if (it != want_.end()) {
      slot_of_[l] = static_cast<std::int32_t>(it->second);
      syms_[it->second].locate = l;
      syms_[it->second].seen = true;
    }
  }
  void trading_state(Locate l, char state) {
    if (Sym* s = sym(l)) {
      s->trading = state == 'T';
    }
  }

  void add(std::uint64_t ts, Locate l, OrderRef ref, char side, Shares n, Price p) {
    book_.add(l, ref, side, n, p);
    after_event(ts, l);
  }
  void execute(std::uint64_t ts, Locate l, OrderRef ref, Shares n, bool printable = true) {
    const auto r = book_.resting(ref);
    book_.execute(l, ref, n);
    if (r && printable && in_grid(ts)) {
      if (Sym* s = sym(r->locate)) {
        const auto v = static_cast<std::int64_t>(std::min(n, r->shares));
        (r->side == Side::buy ? s->sell_vol : s->buy_vol) += v;
      }
    }
    after_event(ts, r ? r->locate : l);
  }
  void cancel(std::uint64_t ts, Locate l, OrderRef ref, Shares n) {
    const auto r = book_.resting(ref);
    book_.cancel(l, ref, n);
    after_event(ts, r ? r->locate : l);
  }
  void remove(std::uint64_t ts, Locate l, OrderRef ref) {
    const auto r = book_.resting(ref);
    book_.remove(l, ref);
    after_event(ts, r ? r->locate : l);
  }
  void replace(std::uint64_t ts, Locate l, OrderRef old_ref, OrderRef new_ref, Shares n, Price p) {
    const auto r = book_.resting(old_ref);
    book_.replace(l, old_ref, new_ref, n, p);
    after_event(ts, r ? r->locate : l);
  }
  void hidden_trade(std::uint64_t ts, Locate l, char side, Shares n) {
    if (!in_grid(ts)) {
      return;
    }
    if (Sym* s = sym(l)) {
      (side == 'S' ? s->hidden_s : s->hidden_b) += n;
    }
  }

  [[nodiscard]] std::size_t rows() const noexcept { return n_rows_; }
  [[nodiscard]] std::uint64_t grid_time(std::uint64_t k) const noexcept { return o_.start_ns + k * o_.grid_ns; }
  [[nodiscard]] const std::vector<Row>& rows_of(std::size_t i) const noexcept { return syms_[i].rows; }
  [[nodiscard]] bool seen(std::size_t i) const noexcept { return syms_[i].seen; }
  [[nodiscard]] const book::FastBook& book() const noexcept { return book_; }
  [[nodiscard]] const RecorderOptions& options() const noexcept { return o_; }

 private:
  struct Sym {
    Locate locate = 0;
    bool seen = false;
    bool trading = false;
    Top bid;
    Top ask;  // missing levels at UINT32_MAX here, written as 0
    std::array<std::int64_t, kDepth> ofi{};
    std::int64_t buy_vol = 0;
    std::int64_t sell_vol = 0;
    std::int64_t hidden_b = 0;
    std::int64_t hidden_s = 0;
    std::int64_t events = 0;
    bool primed = false;  // bid/ask hold the book's current top
    std::vector<Row> rows;
  };

  [[nodiscard]] Sym* sym(Locate l) noexcept {
    const std::int32_t i = slot_of_[l];
    return i < 0 ? nullptr : &syms_[static_cast<std::size_t>(i)];
  }
  [[nodiscard]] bool in_grid(std::uint64_t ts) const noexcept {
    return ts > o_.start_ns && ts <= o_.start_ns + n_rows_ * o_.grid_ns;
  }

  void read_top(Locate l, Top& bid, Top& ask) const {
    bid = Top{};
    ask = Top{};
    ask.px.fill(std::numeric_limits<std::uint32_t>::max());
    const auto b = book_.levels(l, Side::buy, kDepth);
    for (std::size_t m = 0; m < b.size(); ++m) {
      bid.px[m] = b[m].price;
      bid.sz[m] = static_cast<std::int64_t>(b[m].shares);
    }
    const auto a = book_.levels(l, Side::sell, kDepth);
    for (std::size_t m = 0; m < a.size(); ++m) {
      ask.px[m] = a[m].price;
      ask.sz[m] = static_cast<std::int64_t>(a[m].shares);
    }
  }

  void after_event(std::uint64_t ts, Locate l) {
    Sym* s = sym(l);
    if (s == nullptr) {
      return;
    }
    Top bid;
    Top ask;
    read_top(l, bid, ask);
    if (!s->primed) {
      s->bid = Top{};
      s->ask = Top{};
      s->ask.px.fill(std::numeric_limits<std::uint32_t>::max());
      s->primed = true;
    }
    if (in_grid(ts)) {
      ofi_update(s->bid, bid, s->ask, ask, s->ofi);
      ++s->events;
    }
    s->bid = bid;
    s->ask = ask;
  }

  static std::int32_t clamp32(std::int64_t v) noexcept {
    return static_cast<std::int32_t>(std::clamp<std::int64_t>(v, std::numeric_limits<std::int32_t>::min(),
                                                              std::numeric_limits<std::int32_t>::max()));
  }

  void emit() {
    for (Sym& s : syms_) {
      Row r{};
      for (std::size_t m = 0; m < kDepth; ++m) {
        const bool has_ask = s.ask.px[m] != std::numeric_limits<std::uint32_t>::max();
        r[kBidPx + m] = static_cast<std::int32_t>(s.bid.px[m]);
        r[kBidSz + m] = clamp32(s.bid.sz[m]);
        r[kAskPx + m] = has_ask ? static_cast<std::int32_t>(s.ask.px[m]) : 0;
        r[kAskSz + m] = has_ask ? clamp32(s.ask.sz[m]) : 0;
        r[kOfi + m] = clamp32(s.ofi[m]);
      }
      r[kBuyVol] = clamp32(s.buy_vol);
      r[kSellVol] = clamp32(s.sell_vol);
      r[kHiddenB] = clamp32(s.hidden_b);
      r[kHiddenS] = clamp32(s.hidden_s);
      r[kEvents] = clamp32(s.events);
      r[kTrading] = s.trading ? 1 : 0;
      s.rows.push_back(r);
      s.ofi.fill(0);
      s.buy_vol = s.sell_vol = s.hidden_b = s.hidden_s = s.events = 0;
    }
    ++next_;
  }

  RecorderOptions o_;
  book::FastBook book_;
  std::unordered_map<std::string, std::size_t> want_;
  std::array<std::int32_t, 65536> slot_of_{};
  std::vector<Sym> syms_;
  std::uint64_t n_rows_ = 0;
  std::uint64_t next_ = 1;  // next grid index to emit
};

}  // namespace itchbook::lob
