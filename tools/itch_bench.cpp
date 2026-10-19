// itch-bench: measure the order book's cost on a real trading day.
//
// Input is read (and, for .gz files, decompressed) on a background thread and
// handed over through a lock-free SPSC ring, so the measured thread does only
// framing, decoding and book updates. Its time spent waiting for input is
// measured and subtracted.
//
//   itch-bench [--engine fast|reference] [--latency] [--max-messages N] [--pin CPU] FILE|-
//
// Default: throughput, i.e. wall time of the measuring thread minus input
// waits, divided by messages. --latency also times every book update
// individually into a log-linear histogram and prints percentiles. The timer
// costs a few ns per read (more on ARM), so the throughput figure always comes
// from a run without --latency.

#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <string_view>
#include <type_traits>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

#include "itchbook/bench/clock.hpp"
#include "itchbook/bench/histogram.hpp"
#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/reference_book.hpp"
#include "itchbook/io/source.hpp"
#include "itchbook/io/spsc_ring.hpp"
#include "itchbook/itch/stream.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)
using ull = unsigned long long;

template <class BookT, bool kTime>
struct BenchHandler : itch::NullHandler {
  using NullHandler::on;
  BookT book;
  bench::Histogram hist;
  std::uint64_t updates = 0;

  template <class F>
  void timed(F&& f) {
    ++updates;
    if constexpr (kTime) {
      const std::uint64_t t0 = bench::ticks();
      f();
      hist.record(bench::ticks() - t0);
    } else {
      f();
    }
  }
  void on(const itch::AddOrder& m) {
    timed([&] { book.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price()); });
  }
  void on(const itch::AddOrderMpid& m) {
    timed([&] { book.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price()); });
  }
  void on(const itch::OrderExecuted& m) {
    timed([&] { book.execute(m.stock_locate(), m.order_ref(), m.executed_shares()); });
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    timed([&] { book.execute(m.stock_locate(), m.order_ref(), m.executed_shares()); });
  }
  void on(const itch::OrderCancel& m) {
    timed([&] { book.cancel(m.stock_locate(), m.order_ref(), m.cancelled_shares()); });
  }
  void on(const itch::OrderDelete& m) {
    timed([&] { book.remove(m.stock_locate(), m.order_ref()); });
  }
  void on(const itch::OrderReplace& m) {
    timed([&] { book.replace(m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price()); });
  }
};

// Median cost of two back-to-back timer reads: the floor under every latency sample.
std::uint64_t timer_overhead_ticks() {
  bench::Histogram h;
  for (int i = 0; i < 1'000'000; ++i) {
    const std::uint64_t a = bench::ticks();
    const std::uint64_t b = bench::ticks();
    h.record(b - a);
  }
  return h.quantile(0.5);
}

void pin_this_thread(int cpu) {
  if (cpu < 0) {
    return;
  }
#if defined(__linux__)
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(static_cast<std::size_t>(cpu), &set);
  if (pthread_setaffinity_np(pthread_self(), sizeof set, &set) != 0) {
    std::fprintf(stderr, "itch-bench: could not pin to CPU %d\n", cpu);
  }
#else
  std::fprintf(stderr, "itch-bench: --pin is Linux-only; ignored\n");
#endif
}

template <class BookT, bool kTime>
int run(const std::string& path, const itch::StreamLimits& limits, int pin) {
  const bench::TickCalibration cal;
  auto prefetch = std::make_unique<io::PrefetchSource>(io::open_source(path));
  pin_this_thread(pin);  // after the producer thread exists, so it does not inherit the pin
  // Heap-allocate: the handler holds the book and a 60 KB histogram.
  auto h = std::make_unique<BenchHandler<BookT, kTime>>();
  const auto t0 = std::chrono::steady_clock::now();
  const auto st = itch::parse_stream(*prefetch, *h, limits);
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const double waited = static_cast<double>(prefetch->wait_ns()) / 1e9;
  const double busy = wall - waited;
  const double ns_tick = cal.ns_per_tick();
  const auto& s = h->book.stats();

  std::printf("engine\t%s\n", std::is_same_v<BookT, book::FastBook> ? "fast" : "reference");
  std::printf("file\t%s\n", path.c_str());
  std::printf("messages\t%llu\n", static_cast<ull>(st.messages));
  std::printf("book_updates\t%llu\n", static_cast<ull>(h->updates));
  std::printf("orphans\t%llu\n", static_cast<ull>(s.orphans()));
  std::printf("anomalies\t%llu\n", static_cast<ull>(s.anomalies()));
  std::printf("max_live_orders\t%llu\n", static_cast<ull>(s.max_live_orders));
  std::printf("wall_s\t%.3f\n", wall);
  std::printf("input_wait_s\t%.3f\n", waited);
  std::printf("busy_s\t%.3f\n", busy);
  std::printf("ns_per_message\t%.1f\n", busy * 1e9 / static_cast<double>(st.messages == 0 ? 1 : st.messages));
  std::printf("messages_per_s\t%.0f\n", static_cast<double>(st.messages) / busy);
  if constexpr (kTime) {
    const auto ns = [&](std::uint64_t t) { return static_cast<double>(t) * ns_tick; };
    std::printf("timer\t%s\t%.2f ns/tick\toverhead %.1f ns\n", bench::tick_source(), ns_tick,
                ns(timer_overhead_ticks()));
    std::printf("update_mean_ns\t%.1f\n", h->hist.mean() * ns_tick);
    for (const double q : {0.5, 0.9, 0.99, 0.999, 0.9999}) {
      std::printf("update_p%g_ns\t%.0f\n", q * 100, ns(h->hist.quantile(q)));
    }
    std::printf("update_max_ns\t%.0f\n", ns(h->hist.max()));
  }
  if (!st.clean() && !st.stopped_at_limit) {
    std::fprintf(stderr, "itch-bench: WARNING: stream not clean\n");
    return 2;
  }
  return 0;
}

void usage() {
  std::fprintf(stderr,
               "usage: itch-bench [--engine fast|reference] [--latency] [--max-messages N] [--pin CPU] FILE|-\n");
}

}  // namespace

int main(int argc, char** argv) {
  bool fast = true;
  bool latency = false;
  int pin = -1;
  itch::StreamLimits limits;
  std::string path;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    const bool has_value = i + 1 < argc;
    if (a == "--engine" && has_value) {
      const std::string_view e = argv[++i];
      fast = e == "fast";
      if (!fast && e != "reference") {
        usage();
        return 1;
      }
    } else if (a == "--latency") {
      latency = true;
    } else if (a == "--pin" && has_value) {
      const std::string_view v = argv[++i];
      if (std::from_chars(v.data(), v.data() + v.size(), pin).ec != std::errc{}) {
        usage();
        return 1;
      }
    } else if (a == "--max-messages" && has_value) {
      const std::string_view v = argv[++i];
      if (std::from_chars(v.data(), v.data() + v.size(), limits.max_messages).ec != std::errc{}) {
        usage();
        return 1;
      }
    } else if (path.empty() && (a == "-" || !a.starts_with("-"))) {
      path = a;
    } else {
      usage();
      return 1;
    }
  }
  if (path.empty()) {
    usage();
    return 1;
  }
  try {
    if (fast) {
      return latency ? run<book::FastBook, true>(path, limits, pin) : run<book::FastBook, false>(path, limits, pin);
    }
    return latency ? run<book::ReferenceBook, true>(path, limits, pin)
                   : run<book::ReferenceBook, false>(path, limits, pin);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-bench: %s\n", e.what());
    return 1;
  }
}
