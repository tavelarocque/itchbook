#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include "itchbook/queue/models.hpp"
#include "itchbook/queue/study.hpp"

namespace itchbook::queue {
namespace {

constexpr auto kPess = static_cast<std::size_t>(Model::pessimistic);
constexpr auto kOpt = static_cast<std::size_t>(Model::optimistic);
constexpr auto kUni = static_cast<std::size_t>(Model::uniform);

// ---- models -------------------------------------------------------------------

TEST(QueueModels, CancelRulesOnAKnownQueue) {
  const Calibration cal;
  // 100 ahead, 200 behind, 150 cancelled.
  EXPECT_DOUBLE_EQ(after_cancel(Model::pessimistic, 100, 300, 150, cal), 100);  // all from behind
  EXPECT_DOUBLE_EQ(after_cancel(Model::optimistic, 100, 300, 150, cal), 0);     // 100 from ahead, rest behind
  EXPECT_DOUBLE_EQ(after_cancel(Model::uniform, 100, 300, 150, cal), 50);       // a third from ahead
  // 250 cancelled from a queue with only 200 behind: at least 50 came from ahead.
  EXPECT_DOUBLE_EQ(after_cancel(Model::pessimistic, 100, 300, 250, cal), 50);
  // Cancelling everything leaves nothing ahead in every model.
  for (std::size_t k = 0; k < kModels; ++k) {
    EXPECT_DOUBLE_EQ(after_cancel(static_cast<Model>(k), 100, 300, 300, cal), 0) << kModelName[k];
  }
}

// hftbacktest's ProbQueueModel depth update (src/backtest/models/queue.rs),
// transcribed: p is the probability the decrease is behind the order.
double hbt_front(double front, double prev_qty, double new_qty, double p) {
  const double chg = prev_qty - new_qty;
  const double a = front;
  const double b = prev_qty - front;
  return std::min(a - (1 - p) * chg + std::min(b - p * chg, 0.0), new_qty);
}

TEST(QueueModels, ProbabilisticModelsMatchHftbacktest) {
  const Calibration cal;
  struct Case {
    Model m;
    double (*p)(double a, double b);
  };
  const Case cases[] = {
      {Model::uniform, [](double a, double b) { return b / (a + b); }},
      {Model::power2, [](double a, double b) { return b * b / (b * b + a * a); }},
      {Model::power3, [](double a, double b) { return b * b * b / (b * b * b + a * a * a); }},
      {Model::log, [](double a, double b) { return std::log1p(b) / (std::log1p(b) + std::log1p(a)); }},
      {Model::front3, [](double a, double b) { return 1 - std::pow(a / (a + b), 3); }},
  };
  for (const auto& c : cases) {
    for (const double a : {1.0, 40.0, 100.0, 250.0, 399.0}) {
      for (const double cut : {1.0, 10.0, 75.0, 150.0}) {
        const double others = 400;
        const double hbt = hbt_front(a, others, others - cut, c.p(a, others - a));
        // Where hftbacktest would go below zero, this implementation stops at zero.
        EXPECT_NEAR(after_cancel(c.m, a, others, cut, cal), std::max(0.0, hbt), 1e-9)
            << kModelName[static_cast<std::size_t>(c.m)] << " a=" << a << " c=" << cut;
      }
    }
  }
}

TEST(QueueModels, CalibrationInterpolatesAndPinsTheEnds) {
  std::array<double, Calibration::kBins> phi{};
  phi.fill(0.25);
  const Calibration cal(phi);
  EXPECT_DOUBLE_EQ(cal(0.0), 0.0);
  EXPECT_DOUBLE_EQ(cal(1.0), 1.0);
  EXPECT_DOUBLE_EQ(cal(0.5), 0.25);
  EXPECT_NEAR(cal(Calibration::centre(0) / 2), 0.125, 1e-12);
  const Calibration identity;  // unfitted: the uniform model
  for (const double r : {0.0, 0.1, 0.33, 0.5, 0.9, 1.0}) {
    EXPECT_NEAR(identity(r), r, 1e-12);
  }
}

// Same history: optimistic <= uniform <= pessimistic, and every model stays
// inside [0, others]. Holds because each update is monotone in `ahead`.
RC_GTEST_PROP(QueueModelProperties, OptimisticUniformPessimisticStayOrdered, ()) {
  const Calibration cal;
  double others = *rc::gen::inRange(1, 2000);
  std::array<double, kModels> a{};
  a.fill(others);
  const int n = *rc::gen::inRange(1, 200);
  for (int i = 0; i < n; ++i) {
    const int kind = *rc::gen::inRange(0, 3);
    const double q = *rc::gen::inRange(0, 500);
    if (kind == 0) {
      others += q;  // add behind
    } else if (kind == 1 && others > 0) {
      const double c = std::min(q, others);
      for (std::size_t k = 0; k < kModels; ++k) {
        a[k] = after_cancel(static_cast<Model>(k), a[k], others, c, cal);
      }
      others -= c;
    } else if (others > 0) {
      const double v = std::min(q, others);
      for (double& e : a) {
        e = std::max(0.0, e - v);
      }
      others -= v;
    }
    RC_ASSERT(a[kOpt] <= a[kUni] + 1e-9);
    RC_ASSERT(a[kUni] <= a[kPess] + 1e-9);
    for (const double e : a) {
      RC_ASSERT(e >= 0.0);
      RC_ASSERT(e <= others + 1e-9);
    }
  }
}

// ---- the study ------------------------------------------------------------------

QueueStudy active_study(std::uint64_t sample = 1) {
  StudyOptions o;
  o.sample = sample;
  QueueStudy s(o);
  s.set_active(true, 0);
  for (Locate l = 0; l < 8; ++l) {
    s.set_trading(l, true);
  }
  return s;
}

TEST(QueueStudy, HandWorkedQueue) {
  QueueStudy s = active_study();
  s.add(1, 1, 10, 'B', 100, 5000);  // A: starts the level, so it is not tracked (no queue ahead)
  EXPECT_FALSE(s.true_ahead(10));
  s.add(2, 1, 11, 'B', 50, 5000);  // T joins behind 100
  ASSERT_TRUE(s.true_ahead(11));
  EXPECT_EQ(*s.true_ahead(11), 100);
  s.add(3, 1, 12, 'B', 200, 5000);  // C joins behind T
  EXPECT_EQ(*s.true_ahead(11), 100);
  EXPECT_EQ(*s.others(11), 300u);
  s.remove(4, 1, 12);  // C cancelled: it was behind
  EXPECT_EQ(*s.true_ahead(11), 100);
  auto e = *s.estimates(11);
  EXPECT_DOUBLE_EQ(e[kPess], 100);
  EXPECT_DOUBLE_EQ(e[kOpt], 0);
  EXPECT_NEAR(e[kUni], 100 - 200.0 / 3, 1e-9);
  s.execute(5, 1, 10, 100);  // A fully executed: T is at the head
  EXPECT_EQ(*s.true_ahead(11), 0);
  s.execute(6, 1, 11, 50);  // T filled
  // C also joined the best with a queue ahead, so it was tracked too; it ended when deleted.
  ASSERT_EQ(s.records().size(), 2u);
  EXPECT_EQ(s.records()[0].ref, 12u);
  EXPECT_EQ(s.records()[0].end, End::deleted);
  const OrderRecord& r = s.records()[1];
  EXPECT_EQ(r.ref, 11u);
  EXPECT_EQ(r.end, End::filled);
  EXPECT_EQ(r.ahead0, 100u);
  EXPECT_EQ(r.filled, 50u);
  EXPECT_EQ(r.first_fill_ts, 6u);
  EXPECT_EQ(r.truth_at_fill, 0);
  EXPECT_EQ(r.n_eval, 2u);
  // Optimistic and uniform thought the first print reached T: predicted fill at ts 5.
  EXPECT_EQ(r.m[kOpt].pred_fill_ts, 5u);
  EXPECT_EQ(r.m[kUni].pred_fill_ts, 5u);
  EXPECT_EQ(r.m[kPess].pred_fill_ts, 6u);
  // Errors just before each print: pessimistic was exact both times.
  EXPECT_DOUBLE_EQ(r.m[kPess].sum_abs_err, 0);
  EXPECT_DOUBLE_EQ(r.m[kOpt].sum_abs_err, 100);  // 0 vs 100 before the first print
  EXPECT_NEAR(r.m[kUni].sum_err, -200.0 / 3, 1e-9);
}

TEST(QueueStudy, ReplaceLosesPriorityAndEndsTracking) {
  QueueStudy s = active_study();
  s.add(1, 1, 10, 'S', 100, 7000);
  s.add(2, 1, 11, 'S', 40, 7000);     // tracked, 100 ahead
  s.add(3, 1, 12, 'S', 60, 7000);     // tracked, 140 ahead
  s.replace(4, 1, 11, 20, 40, 7000);  // 11 re-enters at the back as 20
  EXPECT_EQ(*s.true_ahead(12), 100);  // 11 was ahead of 12 and left
  EXPECT_EQ(*s.others(12), 140u);     // net level change for 12: none
  ASSERT_EQ(s.records().size(), 1u);
  EXPECT_EQ(s.records()[0].end, End::replaced);
  s.set_active(false, 9);  // the close ends everything still tracked
  ASSERT_EQ(s.records().size(), 2u);
  EXPECT_EQ(s.records()[1].end, End::close);
}

TEST(QueueStudy, OnlyOrdersJoiningTheBestPriceAreTracked) {
  QueueStudy s = active_study();
  s.add(1, 1, 10, 'B', 100, 5000);
  s.add(2, 1, 11, 'B', 100, 4900);  // behind the best price
  s.add(3, 1, 12, 'B', 100, 5100);  // improves the price: nothing ahead
  s.add(4, 1, 13, 'B', 100, 5100);  // joins the new best
  EXPECT_FALSE(s.true_ahead(11));
  EXPECT_FALSE(s.true_ahead(12));
  ASSERT_TRUE(s.true_ahead(13));
  EXPECT_EQ(*s.true_ahead(13), 100);
  QueueStudy closed = active_study();
  closed.set_active(false, 0);
  closed.add(1, 1, 10, 'B', 100, 5000);
  closed.add(2, 1, 11, 'B', 100, 5000);
  EXPECT_FALSE(closed.true_ahead(11));
}

// Truth against the L3 book after every event: the incremental count of
// shares ahead must equal what FastBook finds by walking the queue.
RC_GTEST_PROP(QueueStudyProperties, TruthMatchesTheBookAfterEveryEvent, ()) {
  constexpr Locate kMaxLocate = 2;
  QueueStudy s = active_study();
  OrderRef next_ref = 1;
  std::vector<OrderRef> issued;
  const int n = *rc::gen::inRange(1, 600);
  for (int i = 0; i < n; ++i) {
    const auto ts = static_cast<std::uint64_t>(i + 1);
    const int kind = *rc::gen::inRange(0, 100);
    const auto l = *rc::gen::inRange<Locate>(1, kMaxLocate + 1);
    const auto price = *rc::gen::inRange<Price>(1, 5);
    const auto shares = *rc::gen::inRange<Shares>(1, 300);
    const OrderRef target = !issued.empty() && *rc::gen::inRange(0, 10) < 9
                                ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())]
                                : next_ref + 1000;
    if (kind < 45) {
      s.add(ts, l, next_ref, *rc::gen::element('B', 'S'), shares, price);
      issued.push_back(next_ref++);
    } else if (kind < 65) {
      s.execute(ts, l, target, shares);
    } else if (kind < 75) {
      s.cancel(ts, l, target, shares);
    } else if (kind < 88) {
      s.remove(ts, l, target);
    } else {
      s.replace(ts, l, target, next_ref, *rc::gen::inRange<Shares>(0, 300), price);
      issued.push_back(next_ref++);
    }
    for (const OrderRef r : issued) {
      const auto truth = s.true_ahead(r);
      if (!truth) {
        continue;
      }
      const auto info = s.book().order(r);
      RC_ASSERT(info.has_value());
      RC_ASSERT(*truth == static_cast<std::int64_t>(info->shares_ahead));
      const auto lv = s.book().levels(info->locate, info->side);
      const auto it = std::find_if(lv.begin(), lv.end(), [&](const auto& x) { return x.price == info->price; });
      RC_ASSERT(it != lv.end());
      RC_ASSERT(*s.others(r) == it->shares - info->shares);
      const auto e = *s.estimates(r);
      RC_ASSERT(e[kOpt] <= e[kUni] + 1e-9);
      RC_ASSERT(e[kUni] <= e[kPess] + 1e-9);
      for (const double x : e) {
        RC_ASSERT(x >= 0.0);
        RC_ASSERT(x <= static_cast<double>(*s.others(r)) + 1e-9);
      }
    }
  }
  s.set_active(false, 1'000'000);
  RC_ASSERT(s.tracked_now() == 0u);
  for (const OrderRecord& r : s.records()) {
    RC_ASSERT(r.filled <= r.own0);
    RC_ASSERT(r.n_eval <= r.n_trades);
    for (const ModelResult& m : r.m) {
      RC_ASSERT(m.pred_filled <= static_cast<double>(r.own0) + 1e-9);
    }
  }
}

// The verify mode is the same check run inside the study on real data.
RC_GTEST_PROP(QueueStudyProperties, VerifyModeAgrees, ()) {
  StudyOptions o;
  o.sample = 1;
  o.verify_every = 1;
  QueueStudy s(o);
  s.set_active(true, 0);
  s.set_trading(1, true);
  OrderRef next = 1;
  std::vector<OrderRef> issued;
  const int n = *rc::gen::inRange(1, 400);
  for (int i = 0; i < n; ++i) {
    const auto ts = static_cast<std::uint64_t>(i + 1);
    const auto price = *rc::gen::inRange<Price>(1, 3);
    if (issued.empty() || *rc::gen::inRange(0, 2) == 0) {
      s.add(ts, 1, next, 'B', *rc::gen::inRange<Shares>(1, 100), price);
      issued.push_back(next++);
    } else {
      const OrderRef t = issued[*rc::gen::inRange<std::size_t>(0, issued.size())];
      if (*rc::gen::inRange(0, 2) == 0) {
        s.execute(ts, 1, t, *rc::gen::inRange<Shares>(1, 100));
      } else {
        s.cancel(ts, 1, t, *rc::gen::inRange<Shares>(1, 100));
      }
    }
  }
  RC_ASSERT(s.verify_failures() == 0u);
}

}  // namespace
}  // namespace itchbook::queue
