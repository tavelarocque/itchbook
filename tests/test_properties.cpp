// Property tests (RapidCheck): for arbitrary message sequences and arbitrary
// chunking, the splitter delivers exactly the encoded messages; for arbitrary
// bytes, it never crashes and accounts for every input byte.

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include <cstdint>
#include <random>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"
#include "support/random_stream.hpp"
#include "support/writer.hpp"

namespace itchbook::itch {
namespace {

RC_GTEST_PROP(StreamProperties, ChunkingNeverChangesTheMessages, ()) {
  const auto types = test::known_types();
  const auto picks = *rc::gen::container<std::vector<std::size_t>>(rc::gen::inRange<std::size_t>(0, types.size()));
  const auto seed = *rc::gen::arbitrary<std::uint64_t>();
  const auto chunk = *rc::gen::inRange<std::size_t>(1, 300);
  const auto read_bytes = *rc::gen::inRange<std::size_t>(1, 5000);
  const auto eos_every = *rc::gen::inRange<std::size_t>(0, 20);  // 0 = never insert zero-length frames

  std::mt19937_64 rng(seed);
  std::vector<std::vector<std::uint8_t>> expected;
  std::vector<std::uint8_t> stream;
  std::uint64_t eos = 0;
  for (std::size_t i = 0; i < picks.size(); ++i) {
    expected.push_back(test::random_message(types[picks[i]], rng));
    test::append_frame(stream, expected.back());
    if (eos_every != 0 && i % eos_every == 0) {
      stream.push_back(0);
      stream.push_back(0);
      ++eos;
    }
  }

  io::MemorySource src(stream, chunk);
  test::Recorder rec;
  const auto st = parse_stream(src, rec, {}, read_bytes);
  RC_ASSERT(st.clean());
  RC_ASSERT(st.messages == expected.size());
  RC_ASSERT(st.end_of_session_frames == eos);
  RC_ASSERT(st.bytes == stream.size());
  RC_ASSERT(rec.payloads == expected);
}

RC_GTEST_PROP(StreamProperties, ArbitraryBytesAreFullyAccountedFor, (const std::vector<std::uint8_t>& junk)) {
  const auto chunk = *rc::gen::inRange<std::size_t>(1, 64);
  io::MemorySource src(junk, chunk);
  test::Recorder rec;
  const auto st = parse_stream(src, rec, {}, 128);
  RC_ASSERT(st.bytes + st.truncated_tail_bytes == junk.size());
  RC_ASSERT(st.frames == st.messages + st.unknown_type + st.bad_length + st.end_of_session_frames);
  RC_ASSERT(rec.payloads.size() == st.messages);
  for (const auto& p : rec.payloads) {
    RC_ASSERT(p.size() == kMessageSize[p[0]]);
  }
}

RC_GTEST_PROP(StreamProperties, GzipRoundTrip, ()) {
  const auto types = test::known_types();
  const auto n = *rc::gen::inRange<std::size_t>(0, 400);
  const auto seed = *rc::gen::arbitrary<std::uint64_t>();
  const auto chunk = *rc::gen::inRange<std::size_t>(1, 2000);
  const auto in_buf = *rc::gen::inRange<std::size_t>(1, 2000);
  const auto level = *rc::gen::inRange(0, 10);

  std::mt19937_64 rng(seed);
  std::vector<std::vector<std::uint8_t>> expected;
  std::vector<std::uint8_t> stream;
  for (std::size_t i = 0; i < n; ++i) {
    expected.push_back(test::random_message(types[rng() % types.size()], rng));
    test::append_frame(stream, expected.back());
  }
  const auto gz = test::gzip(stream, level);
  io::GzipSource src(std::make_unique<io::MemorySource>(gz, chunk), in_buf);
  test::Recorder rec;
  const auto st = parse_stream(src, rec, {}, 777);
  RC_ASSERT(st.clean());
  RC_ASSERT(rec.payloads == expected);
}

}  // namespace
}  // namespace itchbook::itch

namespace itchbook::itch {
namespace {

// A handler that also asks for prefetches: delivery must be unchanged, and
// every prefetch must describe a complete frame that is about to be decoded.
struct PrefetchRecorder : test::Recorder {
  std::uint64_t prefetches[2] = {0, 0};
  std::uint64_t bad = 0;
  void prefetch(const std::uint8_t* p, std::size_t len, int stage) {
    ++prefetches[stage];
    // Touch the last byte: ASan flags any frame that runs past the buffer.
    volatile std::uint8_t last = p[len - 1];
    (void)last;
    if (kMessageSize[p[0]] != 0 && kMessageSize[p[0]] != len) {
      ++bad;  // the stream below only contains well-formed frames
    }
  }
};
static_assert(Prefetching<PrefetchRecorder>);
static_assert(!Prefetching<test::Recorder>);

RC_GTEST_PROP(StreamProperties, PrefetchingDoesNotChangeDelivery, ()) {
  const auto types = test::known_types();
  const auto n = *rc::gen::inRange<std::size_t>(0, 3000);
  const auto seed = *rc::gen::arbitrary<std::uint64_t>();
  const auto chunk = *rc::gen::inRange<std::size_t>(1, 400);
  const auto read_bytes = *rc::gen::inRange<std::size_t>(1, 3000);
  std::mt19937_64 rng(seed);
  std::vector<std::vector<std::uint8_t>> expected;
  std::vector<std::uint8_t> stream;
  for (std::size_t i = 0; i < n; ++i) {
    expected.push_back(test::random_message(types[rng() % types.size()], rng));
    test::append_frame(stream, expected.back());
  }
  io::MemorySource src(stream, chunk);
  PrefetchRecorder rec;
  const auto st = parse_stream(src, rec, {}, read_bytes);
  RC_ASSERT(st.clean());
  RC_ASSERT(rec.payloads == expected);
  RC_ASSERT(rec.bad == 0U);
  // Every frame is announced at least once per stage before it is decoded.
  RC_ASSERT(rec.prefetches[0] >= n);
  RC_ASSERT(rec.prefetches[1] >= n);
}

}  // namespace
}  // namespace itchbook::itch
