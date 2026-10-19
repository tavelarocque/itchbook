#pragma once
// A cheap per-event timer for latency histograms.
//
// x86-64: the time-stamp counter (rdtsc), a constant-rate counter on any CPU
// from the last decade, read in a few nanoseconds. ARM64: the generic timer's
// virtual counter (cntvct_el0). On Apple silicon that counter runs at 24 MHz,
// one tick per 41.7 ns, so per-event latencies there are quantised to 41.7 ns
// and the benchmark reports the resolution next to every percentile.
// Elsewhere: steady_clock.
//
// Ticks are converted to nanoseconds with a rate calibrated against
// steady_clock over the whole run, so no CPU frequency is assumed.

#include <chrono>
#include <cstdint>

#if defined(__x86_64__)
#include <x86intrin.h>
#endif

namespace itchbook::bench {

[[nodiscard]] inline std::uint64_t ticks() noexcept {
#if defined(__x86_64__)
  return __rdtsc();
#elif defined(__aarch64__)
  std::uint64_t v;
  asm volatile("mrs %0, cntvct_el0" : "=r"(v));
  return v;
#else
  return static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count());
#endif
}

[[nodiscard]] inline const char* tick_source() noexcept {
#if defined(__x86_64__)
  return "rdtsc";
#elif defined(__aarch64__)
  return "cntvct_el0";
#else
  return "steady_clock";
#endif
}

class TickCalibration {
 public:
  TickCalibration() : t0_(ticks()), c0_(std::chrono::steady_clock::now()) {}
  // Nanoseconds per tick measured from construction until now.
  [[nodiscard]] double ns_per_tick() const {
    const auto dt = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - c0_).count();
    const auto dk = static_cast<double>(ticks() - t0_);
    return dk > 0 ? dt / dk : 1.0;
  }

 private:
  std::uint64_t t0_;
  std::chrono::steady_clock::time_point c0_;
};

}  // namespace itchbook::bench
