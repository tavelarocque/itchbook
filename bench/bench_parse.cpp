// Parser microbenchmarks on an in-memory synthetic day slice whose message
// mix matches the published 01302020 per-type counts. These isolate framing
// and decoding from disk and gzip; full-file numbers come from itch-stats.

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/fields.hpp"
#include "itchbook/itch/stream.hpp"
#include "support/random_stream.hpp"
#include "support/writer.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)

// Message mix of 01302020 (counts per million messages, rounded).
const std::vector<std::pair<char, int>> kMix = {{'A', 436436}, {'D', 425920}, {'U', 86886}, {'E', 19882}, {'X', 11791},
                                                {'I', 9509},   {'F', 4430},   {'P', 4205},  {'L', 512},   {'C', 329},
                                                {'Q', 42},     {'Y', 21},     {'H', 21},    {'R', 21}};

const std::vector<std::uint8_t>& day_slice() {
  static const std::vector<std::uint8_t> s = [] {
    std::mt19937_64 rng(42);
    std::vector<char> bag;
    for (const auto& [t, w] : kMix) {
      bag.insert(bag.end(), static_cast<std::size_t>(w), t);
    }
    std::vector<std::uint8_t> out;
    for (int i = 0; i < 2'000'000; ++i) {
      test::append_frame(out, test::random_message(bag[rng() % bag.size()], rng));
    }
    return out;
  }();
  return s;
}

struct CountOnly : itch::NullHandler {
  using NullHandler::on;
  std::uint64_t n = 0;
  template <class M>
  void on(const M& /*m*/) {
    ++n;
  }
};

// Decodes every field of every message, so nothing can be optimized away.
struct TouchAll : itch::NullHandler {
  using NullHandler::on;
  std::uint64_t acc = 0;
  template <class M>
  void on(const M& m) {
    itch::visit_fields(m, [&](std::string_view, std::uint64_t v) { acc += v; });
  }
};

// Reads what a book update needs from the hot messages.
struct BookFields : itch::NullHandler {
  using NullHandler::on;
  std::uint64_t acc = 0;
  void on(const itch::AddOrder& m) {
    acc += m.order_ref() + m.shares() + m.price() + static_cast<std::uint8_t>(m.side());
  }
  void on(const itch::OrderDelete& m) { acc += m.order_ref(); }
  void on(const itch::OrderReplace& m) { acc += m.original_order_ref() + m.new_order_ref() + m.shares() + m.price(); }
  void on(const itch::OrderExecuted& m) { acc += m.order_ref() + m.executed_shares(); }
  void on(const itch::OrderCancel& m) { acc += m.order_ref() + m.cancelled_shares(); }
};

template <class H>
void BM_Parse(benchmark::State& state) {
  const auto& s = day_slice();
  std::uint64_t msgs = 0;
  for (auto _ : state) {
    io::MemorySource src(s);
    H h;
    const auto st = itch::parse_stream(src, h);
    benchmark::DoNotOptimize(h);
    msgs += st.messages;
  }
  state.SetItemsProcessed(static_cast<std::int64_t>(msgs));
  state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(s.size()));
}

BENCHMARK(BM_Parse<CountOnly>)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Parse<BookFields>)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Parse<TouchAll>)->Unit(benchmark::kMillisecond);

void BM_GzipDecode(benchmark::State& state) {
  static const auto gz = test::gzip(day_slice());
  std::vector<std::uint8_t> out(1U << 20U);
  for (auto _ : state) {
    io::GzipSource src(std::make_unique<io::MemorySource>(gz));
    while (src.read(out.data(), out.size()) != 0) {
    }
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(day_slice().size()));
}
BENCHMARK(BM_GzipDecode)->Unit(benchmark::kMillisecond);

}  // namespace

BENCHMARK_MAIN();
