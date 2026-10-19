// Benchmark support: the log-linear histogram and the prefetching SPSC source.

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "itchbook/bench/histogram.hpp"
#include "itchbook/io/source.hpp"
#include "itchbook/io/spsc_ring.hpp"
#include "itchbook/itch/stream.hpp"
#include "support/random_stream.hpp"

namespace itchbook {
namespace {

TEST(Histogram, BucketsRoundTripAndBoundRelativeError) {
  for (std::uint64_t v : {0ULL, 1ULL, 255ULL, 256ULL, 257ULL, 511ULL, 512ULL, 1000ULL, 123456789ULL, 1ULL << 62U}) {
    const auto i = bench::Histogram::index(v);
    const auto lo = bench::Histogram::lower_bound(i);
    EXPECT_LE(lo, v);
    EXPECT_EQ(bench::Histogram::index(lo), i);
    EXPECT_LE(static_cast<double>(v - lo), static_cast<double>(v) / 128.0 + 1.0) << v;
  }
  EXPECT_LT(bench::Histogram::index(~0ULL), bench::Histogram::kSize);
}

RC_GTEST_PROP(HistogramProperties, QuantilesMatchSortedSamplesWithinOnePercent,
              (const std::vector<std::uint32_t>& raw)) {
  RC_PRE(!raw.empty());
  bench::Histogram h;
  std::vector<std::uint64_t> v(raw.begin(), raw.end());
  for (const auto x : v) {
    h.record(x);
  }
  std::sort(v.begin(), v.end());
  for (const double q : {0.5, 0.9, 0.99, 1.0}) {
    auto rank = static_cast<std::size_t>(q * static_cast<double>(v.size()));
    rank = std::clamp<std::size_t>(rank, 1, v.size());
    const std::uint64_t exact = v[rank - 1];
    const std::uint64_t got = h.quantile(q);
    RC_ASSERT(got <= exact);
    RC_ASSERT(static_cast<double>(exact - got) <= static_cast<double>(exact) / 128.0 + 1.0);
  }
}

RC_GTEST_PROP(PrefetchProperties, DeliversExactlyTheInnerBytes, ()) {
  const auto n = *rc::gen::inRange<std::size_t>(0, 200'000);
  const auto chunk = *rc::gen::inRange<std::size_t>(1, 5000);
  const auto blocks = *rc::gen::inRange<std::size_t>(1, 6);
  const auto block_bytes = *rc::gen::inRange<std::size_t>(1, 3000);
  const auto read_cap = *rc::gen::inRange<std::size_t>(1, 7000);
  std::vector<std::uint8_t> data(n);
  std::mt19937_64 rng(n);
  for (auto& b : data) {
    b = static_cast<std::uint8_t>(rng());
  }
  io::PrefetchSource src(std::make_unique<io::MemorySource>(data, chunk), blocks, block_bytes);
  std::vector<std::uint8_t> out;
  std::vector<std::uint8_t> buf(read_cap);
  for (;;) {
    const std::size_t got = src.read(buf.data(), buf.size());
    if (got == 0) {
      break;
    }
    out.insert(out.end(), buf.begin(), buf.begin() + static_cast<std::ptrdiff_t>(got));
  }
  RC_ASSERT(out == data);
}

TEST(Prefetch, EarlyDestructionWithAFullRingDoesNotHang) {
  std::vector<std::uint8_t> data(1U << 20U, 7);
  for (int i = 0; i < 20; ++i) {
    io::PrefetchSource src(std::make_unique<io::MemorySource>(data, 100), 2, 64);
    std::uint8_t b[10];
    (void)src.read(b, sizeof b);
  }  // destructor must stop the producer
  SUCCEED();
}

TEST(Prefetch, ParsesAGzipStreamIdenticallyToDirectReading) {
  std::mt19937_64 rng(5);
  const auto types = test::known_types();
  std::vector<std::vector<std::uint8_t>> msgs;
  std::vector<std::uint8_t> raw;
  for (int i = 0; i < 20000; ++i) {
    msgs.push_back(test::random_message(types[rng() % types.size()], rng));
    test::append_frame(raw, msgs.back());
  }
  const auto gz = test::gzip(raw);
  io::PrefetchSource src(std::make_unique<io::GzipSource>(std::make_unique<io::MemorySource>(gz, 999)), 4, 4096);
  test::Recorder rec;
  const auto st = itch::parse_stream(src, rec);
  EXPECT_TRUE(st.clean());
  EXPECT_EQ(rec.payloads, msgs);
}

}  // namespace
}  // namespace itchbook
