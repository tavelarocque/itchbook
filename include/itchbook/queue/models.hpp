#pragma once
// L2 queue-position models: estimates of how many shares are ahead of a
// resting order, using only what a market-by-price feed shows.
//
// The observer sees the size of the other orders at its price level (the
// level without its own order, which is what a backtest that simulates the
// order would see) and every trade print at that price. A trade always
// consumes the front of the queue. A decrease without a trade is a
// cancellation, and the models differ only in where they assume it came from.
// Every model has the form
//
//     ahead' = ahead - c * phi(ahead, behind)
//
// for a cancellation of c shares, clamped to what is arithmetically possible:
// at most c shares can leave from ahead, and at most `behind` from behind, so
//
//     max(0, ahead - c) <= ahead' <= min(ahead, (ahead + behind) - c).
//
// The two clamps are the optimistic and pessimistic models themselves.
//
// The probabilistic models are those of hftbacktest (nkaz001/hftbacktest,
// src/backtest/models/queue.rs), which writes p for the probability that a
// decrease is BEHIND the order, so phi = 1 - p here:
//   pessimistic  RiskAdverseQueueModel
//   power2/3     ProbQueueModel + PowerProbQueueFunc(n):  p = B^n / (B^n + A^n)
//   log          ProbQueueModel + LogProbQueueFunc:       p = ln(1+B) / (ln(1+B) + ln(1+A))
//   front3       ProbQueueModel + PowerProbQueueFunc3(3): p = 1 - (A/(A+B))^3, the model
//                hftbacktest's own Level-3 example compares against
// Two differences from hftbacktest, both deliberate: an estimate never goes
// below zero on a cancellation (hftbacktest lets `front` go negative, which
// later reads as a fill), and A and B are both taken after the trades already
// applied (hftbacktest's B includes shares traded since the last depth update).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>

namespace itchbook::queue {

enum class Model : std::uint8_t {
  pessimistic,  // cancellations come from behind first: only trades advance the order
  optimistic,   // cancellations come from ahead first
  uniform,      // a cancelled share is equally likely to be any share in the queue: phi = A/(A+B)
  power2,       // phi = A^2/(A^2+B^2): cancellations concentrate in the larger part of the queue
  power3,       // phi = A^3/(A^3+B^3)
  log,          // phi = log(1+A)/(log(1+A)+log(1+B))
  front3,       // phi = (A/(A+B))^3: cancellations mostly from behind
  calibrated,   // phi = empirical share of cancelled volume ahead, by relative position A/(A+B)
};
inline constexpr std::size_t kModels = 8;
inline constexpr std::array<std::string_view, kModels> kModelName = {
    "pessimistic", "optimistic", "uniform", "power2", "power3", "log", "front3", "calibrated"};

// Empirical phi(r), r = ahead / (ahead + behind), tabulated at the centres of
// equal-width bins over [0, 1] and interpolated linearly, with phi(0) = 0 and
// phi(1) = 1 forced (nothing can be cancelled from an empty side).
class Calibration {
 public:
  static constexpr std::size_t kBins = 20;
  Calibration() {
    for (std::size_t i = 0; i < kBins; ++i) {
      phi_[i] = centre(i);  // uniform until fitted
    }
  }
  explicit Calibration(std::span<const double> phi) {
    for (std::size_t i = 0; i < kBins; ++i) {
      phi_[i] = i < phi.size() ? std::clamp(phi[i], 0.0, 1.0) : centre(i);
    }
    fitted_ = true;
  }
  [[nodiscard]] static constexpr double centre(std::size_t i) noexcept {
    return (static_cast<double>(i) + 0.5) / static_cast<double>(kBins);
  }
  [[nodiscard]] static std::size_t bin(double r) noexcept {
    const auto b = static_cast<std::size_t>(std::max(0.0, r) * static_cast<double>(kBins));
    return std::min(b, kBins - 1);
  }
  [[nodiscard]] double operator()(double r) const noexcept {
    r = std::clamp(r, 0.0, 1.0);
    // Knots: (0,0), (centre(i), phi_i) for each bin, (1,1).
    double x0 = 0.0;
    double y0 = 0.0;
    for (std::size_t i = 0; i <= kBins; ++i) {
      const double x1 = i < kBins ? centre(i) : 1.0;
      const double y1 = i < kBins ? phi_[i] : 1.0;
      if (r <= x1) {
        return x1 > x0 ? y0 + (y1 - y0) * (r - x0) / (x1 - x0) : y1;
      }
      x0 = x1;
      y0 = y1;
    }
    return 1.0;
  }
  [[nodiscard]] bool fitted() const noexcept { return fitted_; }
  [[nodiscard]] const std::array<double, kBins>& table() const noexcept { return phi_; }

 private:
  std::array<double, kBins> phi_{};
  bool fitted_ = false;
};

// Fraction of a cancellation the model attributes to the shares ahead.
[[nodiscard]] inline double phi(Model m, double ahead, double behind, const Calibration& cal) noexcept {
  const double total = ahead + behind;
  if (total <= 0.0) {
    return 0.0;
  }
  switch (m) {
    case Model::pessimistic: return 0.0;
    case Model::optimistic: return 1.0;
    case Model::uniform: return ahead / total;
    case Model::power2: {
      const double a = ahead * ahead;
      const double b = behind * behind;
      return a / (a + b);
    }
    case Model::power3: {
      const double a = ahead * ahead * ahead;
      const double b = behind * behind * behind;
      return a / (a + b);
    }
    case Model::log: {
      const double a = std::log1p(ahead);
      const double b = std::log1p(behind);
      return a + b > 0.0 ? a / (a + b) : 0.0;
    }
    case Model::front3: {
      const double r = ahead / total;
      return r * r * r;
    }
    case Model::calibrated: return cal(ahead / total);
  }
  return 0.0;
}

// A cancellation of c shares from a queue (excluding the order itself) of
// `others` shares, `ahead` of which the model believes are ahead.
[[nodiscard]] inline double after_cancel(Model m, double ahead, double others, double c,
                                         const Calibration& cal) noexcept {
  c = std::min(c, others);
  ahead = std::clamp(ahead, 0.0, others);
  const double est = ahead - c * phi(m, ahead, others - ahead, cal);
  const double lo = std::max(0.0, ahead - c);
  const double hi = std::min(ahead, others - c);
  return std::clamp(est, lo, hi);
}

}  // namespace itchbook::queue
