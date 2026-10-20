#pragma once
// QueueStudy: follows a sample of resting orders through the day and, at every
// trade print at their price, compares each L2 model's estimate of the shares
// ahead with the true number from the L3 book.
//
// Ground truth. Every live order carries its arrival sequence number. Orders
// only ever join a queue at the back, so an order X is ahead of a tracked
// order T exactly when they share a price level and X arrived first. The true
// shares ahead of T therefore change only when an order that arrived earlier
// at T's level loses shares (executed, cancelled, deleted or replaced away),
// and the study keeps that number incrementally, without walking queues. The
// tests and the --verify mode check it against FastBook::order(), which walks
// the queue.
//
// The L2 view. For each tracked order the models see what a market-by-price
// feed shows about the rest of its level: the size of the other orders at the
// price (net change per message) and every trade print there, including
// prints against the tracked order itself, which they cannot tell apart from
// others. They never see the tracked order's own fills.
//
// Population: orders added (A or F) during regular market hours, while the
// symbol is trading, at the current best price on their side when that level
// already has other orders resting, so the order starts with a queue ahead of
// it. A deterministic hash of the order reference selects 1 in `sample`.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/ref_index.hpp"
#include "itchbook/book/types.hpp"
#include "itchbook/queue/models.hpp"

namespace itchbook::queue {

using book::Locate;
using book::OrderRef;
using book::Price;
using book::Shares;
using book::Side;

enum class End : std::uint8_t { filled, cancelled, deleted, replaced, close };
inline constexpr std::array<const char*, 5> kEndName = {"filled", "cancelled", "deleted", "replaced", "close"};

struct ModelResult {
  double sum_abs_err = 0;              // over evaluation points, |estimate - truth| in shares
  double sum_err = 0;                  // signed: positive means the model thinks more is ahead (pessimistic)
  double err_at_fill = 0;              // estimate - truth at the order's first real fill
  std::uint64_t pred_fill_ts = 0;      // first trade print the model would have filled the order on (0 = never)
  std::uint64_t pred_fill_trades = 0;  // trade prints at the level up to and including that one
  double pred_filled = 0;              // shares the model would have filled over the order's life
};

struct OrderRecord {
  OrderRef ref;
  Locate locate;
  Side side;
  Price price;
  std::uint64_t add_ts;
  std::uint64_t end_ts;
  End end;
  Shares own0;
  std::uint64_t ahead0;         // shares ahead at entry (all of the level's other shares)
  std::uint32_t orders_ahead0;  // orders ahead at entry
  Price spread;                 // best ask - best bid at entry, in ITCH price units (1/10000 dollar)
  std::uint64_t first_fill_ts;  // 0 = never filled
  std::uint64_t fill_trades;    // trade prints at the level up to and including the first real fill
  std::uint64_t filled;         // shares actually executed
  std::int64_t truth_at_fill;   // true shares ahead at the first real fill (0 when filled at the head)
  std::uint64_t n_eval;         // evaluation points (trade prints before or at the first real fill)
  std::uint64_t n_trades;       // all trade prints at the level while tracked
  std::array<ModelResult, kModels> m;
};

// One evaluation point of an event-sampled order: the truth and every
// model's estimate just before a trade print at its level.
struct EvalPoint {
  OrderRef ref;
  std::uint64_t ts;
  std::uint64_t others;
  std::int64_t truth;
  std::array<double, kModels> est;
};

// Where cancellations come from: cancelled volume at tracked orders' levels,
// by the tracked order's true relative position r = ahead / others, split into
// volume that was ahead of it and volume behind.
struct CancelLocation {
  std::array<double, Calibration::kBins> ahead{};
  std::array<double, Calibration::kBins> total{};
  std::array<std::uint64_t, Calibration::kBins> events{};
};

// Where the queue around tracked orders goes: shares that left from ahead of
// a tracked order and from behind it, by cause, during regular hours. Summed
// over tracked orders, so a share counts once per tracked order it was ahead
// of or behind.
enum class Cause : std::uint8_t { trade, cancel, del, replace_away, replace_same };
inline constexpr std::size_t kCauses = 5;
inline constexpr std::array<const char*, kCauses> kCauseName = {"trade", "partial_cancel", "delete",
                                                                "replace_to_other_price", "replace_same_price"};
struct QueueFlow {
  std::array<double, kCauses> ahead{};
  std::array<double, kCauses> behind{};
  // Same-price replaces as the L2 feed shows them: the net size change.
  double replace_same_net_visible = 0;  // |net change| summed, for replaces by orders ahead
};

struct StudyOptions {
  std::uint64_t sample = 16;       // track 1 in `sample` eligible orders
  std::uint64_t event_sample = 0;  // of those, record evaluation points for 1 in `event_sample` (0 = none)
  std::uint64_t verify_every = 0;  // check truth against the L3 book every N evaluation points (0 = never)
  Calibration calibration;
  // Sinks. When unset, records and points are kept in memory (for tests).
  std::function<void(const OrderRecord&)> on_record;
  std::function<void(const EvalPoint&)> on_point;
};

class QueueStudy {
 public:
  explicit QueueStudy(StudyOptions o = {}) : o_(std::move(o)) {}

  // Tracking starts only while active (regular hours); the trading state is per symbol.
  void set_active(bool on, std::uint64_t ts) {
    active_ = on;
    if (!on) {
      close_all(ts);
    }
  }
  void set_trading(Locate l, bool trading) { trading_[l] = trading; }

  void add(std::uint64_t ts, Locate locate, OrderRef ref, char side_c, Shares shares, Price price) {
    if ((side_c != 'B' && side_c != 'S') || index_.find(ref) != book::kNoSlot) {
      book_.add(locate, ref, side_c, shares, price);  // the book counts the anomaly
      return;
    }
    const Side side = side_c == 'B' ? Side::buy : Side::sell;
    const std::uint64_t key = level_key(locate, side, price);
    const auto best = book_.best(locate, side);
    const auto other = book_.best(locate, side == Side::buy ? Side::sell : Side::buy);
    const bool eligible =
        active_ && trading_[locate] && best && best->price == price && best->shares > 0 && sampled(ref, o_.sample);
    for_tracked(key, [&](Tracked& t) { t.others += shares; });  // joins behind everyone
    const std::uint32_t s = alloc_order({next_seq_++, price, shares, locate, side, kNone});
    index_.insert(ref, s);
    if (eligible) {
      start_tracking(ts, ref, s, key, *best, other ? other->price : 0);
    }
    book_.add(locate, ref, side_c, shares, price);
  }

  void execute(std::uint64_t ts, Locate locate, OrderRef ref, Shares n) {
    const std::uint32_t s = index_.find(ref);
    if (s != book::kNoSlot) {
      const Ord& x = ords_[s];
      const Shares v = std::min(n, x.shares);
      const std::uint64_t key = level_key(x.locate, x.side, x.price);
      for_tracked(key, [&](Tracked& t) { trade(ts, t, x, t.slot == s, v); });
      reduce(ts, s, ref, v, End::filled);
    }
    book_.execute(locate, ref, n);
  }

  void cancel(std::uint64_t ts, Locate locate, OrderRef ref, Shares n) {
    const std::uint32_t s = index_.find(ref);
    if (s != book::kNoSlot) {
      const Ord& x = ords_[s];
      const Shares c = std::min(n, x.shares);
      const std::uint64_t key = level_key(x.locate, x.side, x.price);
      for_tracked(key, [&](Tracked& t) {
        if (t.slot != s) {
          other_removed(t, x, c, Cause::cancel);
        }
      });
      reduce(ts, s, ref, c, End::cancelled);
    }
    book_.cancel(locate, ref, n);
  }

  void remove(std::uint64_t ts, Locate locate, OrderRef ref) {
    const std::uint32_t s = index_.find(ref);
    if (s != book::kNoSlot) {
      const Ord& x = ords_[s];
      const std::uint64_t key = level_key(x.locate, x.side, x.price);
      for_tracked(key, [&](Tracked& t) {
        if (t.slot != s) {
          other_removed(t, x, x.shares, Cause::del);
        }
      });
      erase(ts, s, ref, End::deleted);
    }
    book_.remove(locate, ref);
  }

  void replace(std::uint64_t ts, Locate locate, OrderRef old_ref, OrderRef new_ref, Shares shares, Price price) {
    const std::uint32_t s = index_.find(old_ref);
    if (s == book::kNoSlot || (new_ref != old_ref && index_.find(new_ref) != book::kNoSlot)) {
      book_.replace(locate, old_ref, new_ref, shares, price);
      return;
    }
    const Ord x = ords_[s];
    const std::uint64_t old_key = level_key(x.locate, x.side, x.price);
    const std::uint64_t new_key = level_key(x.locate, x.side, price);
    if (old_key == new_key) {
      // One message, one level: the L2 feed shows only the net change.
      for_tracked(old_key, [&](Tracked& t) {
        if (t.slot != s) {
          other_replaced_same_level(t, x, shares);
        }
      });
    } else {
      for_tracked(old_key, [&](Tracked& t) {
        if (t.slot != s) {
          other_removed(t, x, x.shares, Cause::replace_away);
        }
      });
      if (shares > 0) {
        for_tracked(new_key, [&](Tracked& t) { t.others += shares; });
      }
    }
    erase(ts, s, old_ref, End::replaced);
    if (shares > 0) {
      const std::uint32_t ns = alloc_order({next_seq_++, price, shares, x.locate, x.side, kNone});
      index_.insert(new_ref, ns);
    }
    book_.replace(locate, old_ref, new_ref, shares, price);
  }

  void close_all(std::uint64_t ts) {
    while (!live_tracked_.empty()) {
      finish(ts, live_tracked_.back(), End::close);
    }
  }

  [[nodiscard]] std::vector<OrderRecord>& records() noexcept { return records_; }
  [[nodiscard]] std::vector<EvalPoint>& eval_points() noexcept { return points_; }
  [[nodiscard]] const CancelLocation& cancel_location() const noexcept { return cancel_loc_; }
  [[nodiscard]] const QueueFlow& queue_flow() const noexcept { return flow_; }
  [[nodiscard]] const book::FastBook& book() const noexcept { return book_; }
  [[nodiscard]] std::uint64_t verified() const noexcept { return verified_; }
  [[nodiscard]] std::uint64_t verify_failures() const noexcept { return verify_failures_; }
  [[nodiscard]] std::size_t tracked_now() const noexcept { return live_tracked_.size(); }

  // Truth for a tracked order, for tests: nullopt when the order is not tracked.
  [[nodiscard]] std::optional<std::int64_t> true_ahead(OrderRef ref) const {
    const std::uint32_t s = index_.find(ref);
    if (s == book::kNoSlot || ords_[s].tracked == kNone) {
      return std::nullopt;
    }
    return tracked_[ords_[s].tracked].truth;
  }
  [[nodiscard]] std::optional<std::array<double, kModels>> estimates(OrderRef ref) const {
    const std::uint32_t s = index_.find(ref);
    if (s == book::kNoSlot || ords_[s].tracked == kNone) {
      return std::nullopt;
    }
    return tracked_[ords_[s].tracked].est;
  }
  [[nodiscard]] std::optional<std::uint64_t> others(OrderRef ref) const {
    const std::uint32_t s = index_.find(ref);
    if (s == book::kNoSlot || ords_[s].tracked == kNone) {
      return std::nullopt;
    }
    return tracked_[ords_[s].tracked].others;
  }

  [[nodiscard]] static bool sampled(OrderRef ref, std::uint64_t every) noexcept {
    return every <= 1 || book::splitmix64(ref) % every == 0;
  }

 private:
  static constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

  struct Ord {
    std::uint64_t seq;
    Price price;
    Shares shares;
    Locate locate;
    Side side;
    std::uint32_t tracked;  // index into tracked_, or kNone
  };

  struct Tracked {
    OrderRecord rec;
    std::uint32_t slot;  // the order's slot in ords_
    std::uint64_t seq;
    std::uint64_t key;
    std::uint64_t others;        // shares of the other orders at the level
    std::int64_t truth;          // true shares ahead
    Shares own;                  // the order's remaining shares
    std::uint32_t live_pos = 0;  // position in live_tracked_
    bool filled_yet = false;
    bool events = false;  // record evaluation points
    std::array<double, kModels> est{};
    std::array<double, kModels> pred_own{};  // shares the model believes are still unfilled
  };

  [[nodiscard]] static std::uint64_t level_key(Locate l, Side s, Price p) noexcept {
    return (std::uint64_t{l} << 33U) | (std::uint64_t{static_cast<std::uint8_t>(s)} << 32U) | p;
  }

  template <class F>
  void for_tracked(std::uint64_t key, F&& f) {
    const auto it = at_level_.find(key);
    if (it == at_level_.end()) {
      return;
    }
    for (const std::uint32_t ti : it->second) {
      f(tracked_[ti]);
    }
  }

  std::uint32_t alloc_order(const Ord& o) {
    if (!free_ords_.empty()) {
      const std::uint32_t s = free_ords_.back();
      free_ords_.pop_back();
      ords_[s] = o;
      return s;
    }
    ords_.push_back(o);
    return static_cast<std::uint32_t>(ords_.size() - 1);
  }

  void start_tracking(std::uint64_t ts, OrderRef ref, std::uint32_t s, std::uint64_t key, const book::LevelInfo& best,
                      Price other_best) {
    const Ord& o = ords_[s];
    std::uint32_t ti;
    if (!free_tracked_.empty()) {
      ti = free_tracked_.back();
      free_tracked_.pop_back();
    } else {
      tracked_.emplace_back();
      ti = static_cast<std::uint32_t>(tracked_.size() - 1);
    }
    Tracked& t = tracked_[ti];
    t = Tracked{};
    t.slot = s;
    t.seq = o.seq;
    t.key = key;
    t.others = best.shares;
    t.truth = static_cast<std::int64_t>(best.shares);
    t.own = o.shares;
    t.events = o_.event_sample > 0 && sampled(ref ^ 0x5bd1e995ULL, o_.event_sample);
    t.est.fill(static_cast<double>(best.shares));
    t.pred_own.fill(static_cast<double>(o.shares));
    OrderRecord& r = t.rec;
    r = OrderRecord{};
    r.ref = ref;
    r.locate = o.locate;
    r.side = o.side;
    r.price = o.price;
    r.add_ts = ts;
    r.own0 = o.shares;
    r.ahead0 = best.shares;
    r.orders_ahead0 = best.orders;
    r.spread = other_best == 0 ? 0 : (o.side == Side::buy ? other_best - o.price : o.price - other_best);
    ords_[s].tracked = ti;
    at_level_[key].push_back(ti);
    t.live_pos = static_cast<std::uint32_t>(live_tracked_.size());
    live_tracked_.push_back(ti);
  }

  // A trade print of v shares at t's level, executing order x.
  void trade(std::uint64_t ts, Tracked& t, const Ord& x, bool own, Shares v) {
    ++t.rec.n_trades;
    if (!t.filled_yet) {
      evaluate(ts, t);
    }
    // Models: the print consumes the front of the queue they believe in.
    for (std::size_t k = 0; k < kModels; ++k) {
      const double before = t.est[k];
      const double fill = std::max(0.0, static_cast<double>(v) - before);
      t.est[k] = std::max(0.0, before - static_cast<double>(v));
      if (fill > 0.0 && t.pred_own[k] > 0.0) {
        ModelResult& m = t.rec.m[k];
        if (m.pred_fill_ts == 0) {
          m.pred_fill_ts = ts;
          m.pred_fill_trades = t.rec.n_trades;
        }
        const double f = std::min(fill, t.pred_own[k]);
        m.pred_filled += f;
        t.pred_own[k] -= f;
      }
    }
    if (own) {
      if (!t.filled_yet) {
        t.filled_yet = true;
        t.rec.first_fill_ts = ts;
        t.rec.fill_trades = t.rec.n_trades;
        t.rec.truth_at_fill = t.truth;
        for (std::size_t k = 0; k < kModels; ++k) {
          // est was already advanced by this print; the error that matters is
          // the one evaluate() recorded just before it, kept in last_err_.
          t.rec.m[k].err_at_fill = last_err_[k];
        }
      }
      t.rec.filled += v;
      t.own -= v;
    } else {
      flow(t, x, Cause::trade, v);
      if (x.seq < t.seq) {
        t.truth -= v;
      }
      t.others -= v;
      clamp(t);
    }
  }

  // Another order x at t's level lost c shares without trading.
  void flow(const Tracked& t, const Ord& x, Cause c, double shares) {
    if (active_) {
      (x.seq < t.seq ? flow_.ahead : flow_.behind)[static_cast<std::size_t>(c)] += shares;
    }
  }

  void other_removed(Tracked& t, const Ord& x, Shares c, Cause cause) {
    flow(t, x, cause, c);
    if (active_) {
      const double r = t.others > 0 ? static_cast<double>(t.truth) / static_cast<double>(t.others) : 0.0;
      const std::size_t b = Calibration::bin(r);
      const bool ahead = x.seq < t.seq;
      cancel_loc_.total[b] += c;
      cancel_loc_.ahead[b] += ahead ? c : 0;
      ++cancel_loc_.events[b];
    }
    if (x.seq < t.seq) {
      t.truth -= c;
    }
    l2_cancel_at(t, c);
    t.others -= c;
    clamp(t);
  }

  void other_replaced_same_level(Tracked& t, const Ord& x, Shares new_shares) {
    // Truth: the old order leaves (from wherever it was); the new one joins at the back.
    flow(t, x, Cause::replace_same, x.shares);
    if (x.seq < t.seq) {
      t.truth -= x.shares;
      if (active_) {
        flow_.replace_same_net_visible += std::abs(static_cast<double>(new_shares) - static_cast<double>(x.shares));
      }
    }
    if (new_shares < x.shares) {
      l2_cancel_at(t, x.shares - new_shares);
      t.others -= x.shares - new_shares;
    } else {
      t.others += new_shares - x.shares;
    }
    clamp(t);
  }

  void l2_cancel_at(Tracked& t, Shares c) {
    const auto others = static_cast<double>(t.others);
    for (std::size_t k = 0; k < kModels; ++k) {
      t.est[k] = after_cancel(static_cast<Model>(k), t.est[k], others, static_cast<double>(c), o_.calibration);
    }
  }

  static void clamp(Tracked& t) {
    const auto others = static_cast<double>(t.others);
    for (double& e : t.est) {
      e = std::clamp(e, 0.0, others);
    }
  }

  void evaluate(std::uint64_t ts, Tracked& t) {
    ++t.rec.n_eval;
    const auto truth = static_cast<double>(t.truth);
    for (std::size_t k = 0; k < kModels; ++k) {
      const double e = t.est[k] - truth;
      last_err_[k] = e;
      t.rec.m[k].sum_abs_err += std::abs(e);
      t.rec.m[k].sum_err += e;
    }
    if (t.events) {
      const EvalPoint p{t.rec.ref, ts, t.others, t.truth, t.est};
      if (o_.on_point) {
        o_.on_point(p);
      } else {
        points_.push_back(p);
      }
    }
    if (o_.verify_every > 0 && ++eval_count_ % o_.verify_every == 0) {
      const auto info = book_.order(t.rec.ref);
      ++verified_;
      if (!info || static_cast<std::int64_t>(info->shares_ahead) != t.truth || info->shares != t.own) {
        ++verify_failures_;
      }
    }
  }

  // Order in slot s lost n shares; it ends (with reason `why`) at zero.
  void reduce(std::uint64_t ts, std::uint32_t s, OrderRef ref, Shares n, End why) {
    Ord& x = ords_[s];
    if (x.tracked != kNone && why == End::cancelled) {
      // The owner cancelled part of the order: nothing the models believe is
      // still unfilled can exceed what is really left.
      Tracked& t = tracked_[x.tracked];
      t.own -= n;
      for (double& p : t.pred_own) {
        p = std::min(p, static_cast<double>(t.own));
      }
    }
    if (n >= x.shares) {
      erase(ts, s, ref, why);
    } else {
      x.shares -= n;
    }
  }

  void erase(std::uint64_t ts, std::uint32_t s, OrderRef ref, End why) {
    if (ords_[s].tracked != kNone) {
      finish(ts, ords_[s].tracked, why);
    }
    index_.erase(ref);
    free_ords_.push_back(s);
  }

  void finish(std::uint64_t ts, std::uint32_t ti, End why) {
    Tracked& t = tracked_[ti];
    t.rec.end_ts = ts;
    t.rec.end = why;
    if (o_.on_record) {
      o_.on_record(t.rec);
    } else {
      records_.push_back(t.rec);
    }
    auto& v = at_level_[t.key];
    std::erase(v, ti);
    if (v.empty()) {
      at_level_.erase(t.key);
    }
    const std::uint32_t moved = live_tracked_.back();
    live_tracked_[t.live_pos] = moved;
    tracked_[moved].live_pos = t.live_pos;
    live_tracked_.pop_back();
    ords_[t.slot].tracked = kNone;
    free_tracked_.push_back(ti);
  }

  StudyOptions o_;
  book::FastBook book_;
  book::PagedRefIndex index_;
  std::vector<Ord> ords_;
  std::vector<std::uint32_t> free_ords_;
  std::vector<Tracked> tracked_;
  std::vector<std::uint32_t> free_tracked_;
  std::vector<std::uint32_t> live_tracked_;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> at_level_;
  std::array<bool, 65536> trading_{};
  bool active_ = false;
  std::uint64_t next_seq_ = 0;
  std::array<double, kModels> last_err_{};
  std::vector<OrderRecord> records_;
  std::vector<EvalPoint> points_;
  CancelLocation cancel_loc_;
  QueueFlow flow_;
  std::uint64_t eval_count_ = 0;
  std::uint64_t verified_ = 0;
  std::uint64_t verify_failures_ = 0;
};

}  // namespace itchbook::queue
