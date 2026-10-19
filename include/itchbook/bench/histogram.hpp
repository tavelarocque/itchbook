#pragma once
// Log-linear latency histogram in the style of HdrHistogram: exact counts for
// values below 256, then 128 sub-buckets per power of two, so any recorded
// value is reported within 1/128 (0.8%) of itself. Fixed size (about 60 KB),
// no allocation while recording, values up to 2^63.

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>

namespace itchbook::bench {

class Histogram {
 public:
  static constexpr unsigned kSubBits = 7;                             // 128 sub-buckets per power of two
  static constexpr std::size_t kLinear = std::size_t{2} << kSubBits;  // 256 exact buckets
  static constexpr std::size_t kSize = kLinear + (64 - kSubBits - 1) * (std::size_t{1} << kSubBits);

  void record(std::uint64_t v) noexcept {
    ++counts_[index(v)];
    ++total_;
    max_ = v > max_ ? v : max_;
    sum_ += v;
  }

  [[nodiscard]] static std::size_t index(std::uint64_t v) noexcept {
    if (v < kLinear) {
      return static_cast<std::size_t>(v);
    }
    const unsigned shift = static_cast<unsigned>(std::bit_width(v)) - (kSubBits + 1);  // >= 1
    const std::uint64_t sub = v >> shift;                                              // in [128, 256)
    return kLinear + (shift - 1) * (std::size_t{1} << kSubBits) + static_cast<std::size_t>(sub - (1U << kSubBits));
  }
  // Smallest value that maps to bucket i.
  [[nodiscard]] static std::uint64_t lower_bound(std::size_t i) noexcept {
    if (i < kLinear) {
      return i;
    }
    const std::size_t k = i - kLinear;
    const unsigned shift = static_cast<unsigned>(k >> kSubBits) + 1;
    const std::uint64_t sub = (k & ((std::size_t{1} << kSubBits) - 1)) + (1U << kSubBits);
    return sub << shift;
  }

  // Value at quantile q in [0, 1]: the lower bound of the bucket holding the
  // ceil(q * n)-th smallest sample (so p100 is the bucket of the maximum).
  [[nodiscard]] std::uint64_t quantile(double q) const noexcept {
    if (total_ == 0) {
      return 0;
    }
    auto rank = static_cast<std::uint64_t>(q * static_cast<double>(total_));
    rank = rank == 0 ? 1 : (rank > total_ ? total_ : rank);
    std::uint64_t seen = 0;
    for (std::size_t i = 0; i < kSize; ++i) {
      seen += counts_[i];
      if (seen >= rank) {
        return lower_bound(i);
      }
    }
    return max_;
  }
  [[nodiscard]] std::uint64_t count() const noexcept { return total_; }
  [[nodiscard]] std::uint64_t max() const noexcept { return max_; }
  [[nodiscard]] double mean() const noexcept {
    return total_ == 0 ? 0.0 : static_cast<double>(sum_) / static_cast<double>(total_);
  }

 private:
  std::array<std::uint64_t, kSize> counts_{};
  std::uint64_t total_ = 0;
  std::uint64_t max_ = 0;
  std::uint64_t sum_ = 0;  // fine for latencies: 10^9 samples * 10^9 ticks still fits
};

}  // namespace itchbook::bench
