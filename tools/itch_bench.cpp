// itch-bench: measure the order book's cost on a real trading day.
//
// Input is read (and, for .gz files, decompressed) on a background thread and
// handed over through a lock-free SPSC ring, so the measured thread does only
// framing, decoding and book updates. Its time spent waiting for input is
// measured and subtracted.
//
//   itch-bench [--engine fast|reference] [--latency] [--prefetch | --prefetch-stages MASK]
//   [--inline-input] [--max-messages N] [--pin CPU] FILE|-
//
// Default: throughput, i.e. wall time of the measuring thread minus input
// waits, divided by messages. --latency also times every book update
// individually into a log-linear histogram and prints percentiles. The timer
// costs a few ns per read (more on ARM), so the throughput figure always comes
// from a run without --latency.
//
// Software prefetching (FastBook only). Off by default: on the final engine
// (paged index, level pool) it measured 7-18 ns per message slower than plain
// decode-and-apply on a Xeon 8488C, because the lookahead pass costs more
// than the misses it hides. Kept for measurement and for the hash index:
//   default             no lookahead: decode and apply each message in turn
//   --prefetch          lookahead cursors with both stages (index entry 16
//                       frames ahead, then the order slot or best level 8 ahead)
//   --prefetch-stages M lookahead cursors; bit 0 enables stage 0, bit 1
//                       stage 1, so 0 measures the cost of the lookahead alone
//   --no-prefetch       the default, accepted for older scripts

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

template <class BookT, bool kTime, bool kScout = true>
struct BenchHandler : itch::NullHandler {
  using NullHandler::on;
  BookT book;
  bench::Histogram hist;
  std::uint64_t updates = 0;
  unsigned prefetch_stages = 3;  // bit i enables prefetch stage i

  void prefetch(const std::uint8_t* p, std::size_t len, int stage) const noexcept
    requires(kScout && requires(const BookT& b, const std::uint8_t* q, std::size_t n) { b.prefetch(q, n, 0); })
  {
    if (((prefetch_stages >> static_cast<unsigned>(stage)) & 1U) != 0) {
      book.prefetch(p, len, stage);
    }
  }

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

struct Options {
  int pin = -1;
  unsigned prefetch_stages = 3;
  bool scout = false;         // run the parser's lookahead cursors at all
  bool inline_input = false;  // read input on the measured thread, as CppTrader does
};

template <class BookT, bool kTime, bool kScout = true>
int run(const std::string& path, const itch::StreamLimits& limits, const Options& o) {
  const bench::TickCalibration cal;
  std::unique_ptr<io::ByteSource> input;
  io::PrefetchSource* prefetcher = nullptr;
  if (o.inline_input) {
    input = io::open_source(path);
  } else {
    auto p = std::make_unique<io::PrefetchSource>(io::open_source(path));
    prefetcher = p.get();
    input = std::move(p);
  }
  pin_this_thread(o.pin);  // after the producer thread exists, so it does not inherit the pin
  // Heap-allocate: the handler holds the book and a 60 KB histogram.
  auto h = std::make_unique<BenchHandler<BookT, kTime, kScout>>();
  h->prefetch_stages = o.prefetch_stages;
  const auto t0 = std::chrono::steady_clock::now();
  const auto st = itch::parse_stream(*input, *h, limits);
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const double waited = prefetcher != nullptr ? static_cast<double>(prefetcher->wait_ns()) / 1e9 : 0.0;
  const double busy = wall - waited;
  const double ns_tick = cal.ns_per_tick();
  const auto& s = h->book.stats();

  std::printf("engine\t%s\n", std::is_same_v<BookT, book::FastBook> ? "fast" : "reference");
  std::printf("file\t%s\n", path.c_str());
  std::printf("prefetch_stages\t%u\n", itch::Prefetching<BenchHandler<BookT, kTime, kScout>> ? o.prefetch_stages : 0U);
  std::printf("lookahead_scan\t%s\n", itch::Prefetching<BenchHandler<BookT, kTime, kScout>> ? "on" : "off");
  std::printf("input\t%s\n", o.inline_input ? "inline" : "background thread");
  std::printf("messages\t%llu\n", static_cast<ull>(st.messages));
  std::printf("book_updates\t%llu\n", static_cast<ull>(h->updates));
  std::printf("orphans\t%llu\n", static_cast<ull>(s.orphans()));
  std::printf("anomalies\t%llu\n", static_cast<ull>(s.anomalies()));
  std::printf("max_live_orders\t%llu\n", static_cast<ull>(s.max_live_orders));
  if constexpr (std::is_same_v<BookT, book::FastBook>) {
    std::printf("order_index\t%s\n", BookT::Index::name());
    if constexpr (requires { h->book.index().peak_memory_bytes(); }) {
      std::printf("order_index_peak_mb\t%.1f\n", static_cast<double>(h->book.index().peak_memory_bytes()) / 1e6);
    } else {
      std::printf("order_index_mb\t%.1f\n", static_cast<double>(h->book.index().memory_bytes()) / 1e6);
    }
  }
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
               "usage: itch-bench [--engine fast|reference] [--latency] [--prefetch | --prefetch-stages MASK] "
               "[--inline-input] "
               "[--max-messages N] [--pin CPU] "
               "FILE|-\n");
}

}  // namespace

int main(int argc, char** argv) {
  bool fast = true;
  bool latency = false;
  Options opt;
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
    } else if (a == "--no-prefetch") {
      opt.scout = false;
    } else if (a == "--prefetch") {
      opt.scout = true;
      opt.prefetch_stages = 3;
    } else if (a == "--prefetch-stages" && has_value) {
      opt.scout = true;
      const std::string_view v = argv[++i];
      if (std::from_chars(v.data(), v.data() + v.size(), opt.prefetch_stages).ec != std::errc{}) {
        usage();
        return 1;
      }
    } else if (a == "--inline-input") {
      opt.inline_input = true;
    } else if (a == "--latency") {
      latency = true;
    } else if (a == "--pin" && has_value) {
      const std::string_view v = argv[++i];
      if (std::from_chars(v.data(), v.data() + v.size(), opt.pin).ec != std::errc{}) {
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
      if (!opt.scout) {
        return latency ? run<book::FastBook, true, false>(path, limits, opt)
                       : run<book::FastBook, false, false>(path, limits, opt);
      }
      return latency ? run<book::FastBook, true>(path, limits, opt) : run<book::FastBook, false>(path, limits, opt);
    }
    return latency ? run<book::ReferenceBook, true>(path, limits, opt)
                   : run<book::ReferenceBook, false>(path, limits, opt);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-bench: %s\n", e.what());
    return 1;
  }
}
