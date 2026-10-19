// Frame splitting and input sources: every message is delivered exactly once
// and byte-for-byte intact, however the input is chunked, compressed or cut.

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <random>
#include <string>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"
#include "support/random_stream.hpp"
#include "support/writer.hpp"

namespace itchbook::itch {
namespace {

using test::append_frame;
using test::Recorder;

struct Fixture {
  std::vector<std::vector<std::uint8_t>> messages;
  std::vector<std::uint8_t> stream;
};

Fixture make_fixture(std::size_t n, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  const auto types = test::known_types();
  Fixture f;
  for (std::size_t i = 0; i < n; ++i) {
    f.messages.push_back(test::random_message(types[rng() % types.size()], rng));
    append_frame(f.stream, f.messages.back());
  }
  return f;
}

StreamStats parse_bytes(std::span<const std::uint8_t> bytes, Recorder& rec, std::size_t chunk,
                        std::size_t read_bytes = std::size_t{1} << 20U, StreamLimits limits = {}) {
  io::MemorySource src(bytes, chunk);
  return parse_stream(src, rec, limits, read_bytes);
}

TEST(Stream, EveryChunkSizeYieldsTheSameMessages) {
  const auto f = make_fixture(2000, 1);
  for (const std::size_t chunk : std::initializer_list<std::size_t>{1, 2, 3, 7, 11, 36, 50, 64, 4096, 1U << 20U}) {
    Recorder rec;
    const auto st = parse_bytes(f.stream, rec, chunk);
    EXPECT_TRUE(st.clean()) << chunk;
    EXPECT_EQ(st.messages, f.messages.size()) << chunk;
    EXPECT_EQ(st.bytes, f.stream.size()) << chunk;
    EXPECT_EQ(rec.payloads, f.messages) << chunk;
  }
}

TEST(Stream, TinyReadBufferForcesCompactionEveryFrame) {
  const auto f = make_fixture(500, 2);
  for (const std::size_t read_bytes : std::initializer_list<std::size_t>{1, 2, 37, 100}) {
    Recorder rec;
    const auto st = parse_bytes(f.stream, rec, SIZE_MAX, read_bytes);
    EXPECT_TRUE(st.clean());
    EXPECT_EQ(rec.payloads, f.messages) << read_bytes;
  }
}

TEST(Stream, MaximalFrameFitsInBuffer) {
  // A 65535-byte frame of an unknown type: must be skipped, not overflow.
  std::vector<std::uint8_t> payload(0xFFFF, 0);
  payload[0] = '~';
  std::vector<std::uint8_t> s;
  append_frame(s, payload);
  const auto f = make_fixture(10, 3);
  s.insert(s.end(), f.stream.begin(), f.stream.end());
  Recorder rec;
  const auto st = parse_bytes(s, rec, 1000, 1);
  EXPECT_EQ(st.unknown_type, 1U);
  EXPECT_EQ(st.messages, 10U);
  EXPECT_EQ(rec.payloads, f.messages);
}

TEST(Stream, ZeroLengthFramesAreCountedAndSkipped) {
  auto f = make_fixture(5, 4);
  std::vector<std::uint8_t> s(f.stream.begin(), f.stream.end());
  s.insert(s.begin() + 2 + static_cast<std::ptrdiff_t>(f.messages[0].size()), {0, 0});
  s.push_back(0);
  s.push_back(0);
  Recorder rec;
  const auto st = parse_bytes(s, rec, 3);
  EXPECT_EQ(st.end_of_session_frames, 2U);
  EXPECT_EQ(st.frames, 7U);
  EXPECT_EQ(st.messages, 5U);
  EXPECT_TRUE(st.clean());
  EXPECT_EQ(rec.payloads, f.messages);
}

TEST(Stream, UnknownTypeAndBadLengthAreCountedAndParsingContinues) {
  const auto f = make_fixture(3, 5);
  std::vector<std::uint8_t> s;
  append_frame(s, f.messages[0]);
  append_frame(s, std::vector<std::uint8_t>{'Z', 1, 2, 3});  // unknown type
  append_frame(s, f.messages[1]);
  std::vector<std::uint8_t> short_add(35, 0);  // 'A' needs 36
  short_add[0] = 'A';
  append_frame(s, short_add);
  append_frame(s, f.messages[2]);
  Recorder rec;
  const auto st = parse_bytes(s, rec, 5);
  EXPECT_EQ(st.unknown_type, 1U);
  EXPECT_EQ(st.bad_length, 1U);
  EXPECT_EQ(st.messages, 3U);
  EXPECT_FALSE(st.clean());
  EXPECT_EQ(rec.payloads, f.messages);
}

TEST(Stream, TruncatedFinalFrameIsReportedNotDelivered) {
  const auto f = make_fixture(20, 6);
  for (std::size_t cut = 1; cut <= 2 + f.messages.back().size() - 1; ++cut) {
    std::vector<std::uint8_t> s(f.stream.begin(), f.stream.end() - static_cast<std::ptrdiff_t>(cut));
    Recorder rec;
    const auto st = parse_bytes(s, rec, 7);
    EXPECT_EQ(st.messages, 19U) << cut;
    EXPECT_EQ(st.truncated_tail_bytes, 2 + f.messages.back().size() - cut) << cut;
    EXPECT_EQ(st.bytes + st.truncated_tail_bytes, s.size());
  }
}

TEST(Stream, MaxMessagesStopsExactly) {
  const auto f = make_fixture(100, 7);
  Recorder rec;
  StreamLimits lim;
  lim.max_messages = 37;
  const auto st = parse_bytes(f.stream, rec, 13, std::size_t{1} << 20U, lim);
  EXPECT_TRUE(st.stopped_at_limit);
  EXPECT_EQ(st.messages, 37U);
  ASSERT_EQ(rec.payloads.size(), 37U);
  EXPECT_EQ(rec.payloads.back(), f.messages[36]);
}

// ---- gzip ------------------------------------------------------------------

StreamStats parse_gzip(std::span<const std::uint8_t> gz, Recorder& rec, std::size_t chunk, std::size_t in_buf) {
  io::GzipSource src(std::make_unique<io::MemorySource>(gz, chunk), in_buf);
  return parse_stream(src, rec);
}

TEST(Gzip, SingleMember) {
  const auto f = make_fixture(3000, 8);
  const auto gz = test::gzip(f.stream);
  for (const std::size_t chunk : std::initializer_list<std::size_t>{1, 17, 4096, 1U << 20U}) {
    Recorder rec;
    const auto st = parse_gzip(gz, rec, chunk, 64);
    EXPECT_TRUE(st.clean()) << chunk;
    EXPECT_EQ(rec.payloads, f.messages) << chunk;
  }
}

TEST(Gzip, MultiMemberIsConcatenated) {
  const auto a = make_fixture(700, 9);
  const auto b = make_fixture(900, 10);
  auto gz = test::gzip(a.stream);
  const auto gz2 = test::gzip(b.stream, 1);
  gz.insert(gz.end(), gz2.begin(), gz2.end());
  Recorder rec;
  const auto st = parse_gzip(gz, rec, 333, 100);
  EXPECT_TRUE(st.clean());
  auto expected = a.messages;
  expected.insert(expected.end(), b.messages.begin(), b.messages.end());
  EXPECT_EQ(rec.payloads, expected);
}

TEST(Gzip, TruncatedFileIsFlaggedAndPrefixStillParses) {
  const auto f = make_fixture(5000, 11);
  auto gz = test::gzip(f.stream);
  gz.resize(gz.size() / 2);
  Recorder rec;
  const auto st = parse_gzip(gz, rec, 4096, 4096);
  EXPECT_TRUE(st.source_truncated);
  EXPECT_FALSE(st.clean());
  ASSERT_GT(rec.payloads.size(), 100U);
  ASSERT_LT(rec.payloads.size(), f.messages.size());
  for (std::size_t i = 0; i < rec.payloads.size(); ++i) {
    ASSERT_EQ(rec.payloads[i], f.messages[i]) << i;
  }
}

TEST(Gzip, CorruptDataThrows) {
  const auto f = make_fixture(2000, 12);
  auto gz = test::gzip(f.stream);
  for (std::size_t i = gz.size() / 3; i < gz.size() / 3 + 64; ++i) {
    gz[i] ^= 0x5A;
  }
  Recorder rec;
  EXPECT_THROW(parse_gzip(gz, rec, 4096, 4096), std::runtime_error);
}

// ---- open_source sniffs gzip magic -------------------------------------------

class TempFile {
 public:
  explicit TempFile(const std::vector<std::uint8_t>& bytes) {
    path_ =
        std::filesystem::temp_directory_path() / ("itchbook_test_" + std::to_string(std::random_device{}()) + ".bin");
    std::ofstream(path_, std::ios::binary)
        .write(reinterpret_cast<const char*>(bytes.data()),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
               static_cast<std::streamsize>(bytes.size()));
  }
  ~TempFile() { std::filesystem::remove(path_); }
  TempFile(const TempFile&) = delete;
  TempFile& operator=(const TempFile&) = delete;
  TempFile(TempFile&&) = delete;
  TempFile& operator=(TempFile&&) = delete;
  [[nodiscard]] std::string path() const { return path_.string(); }

 private:
  std::filesystem::path path_;
};

TEST(OpenSource, RawAndGzipFilesGiveIdenticalMessages) {
  const auto f = make_fixture(1500, 13);
  const TempFile raw(f.stream);
  const TempFile gz(test::gzip(f.stream));
  for (const auto* tf : {&raw, &gz}) {
    auto src = io::open_source(tf->path());
    Recorder rec;
    const auto st = parse_stream(*src, rec);
    EXPECT_TRUE(st.clean());
    EXPECT_EQ(rec.payloads, f.messages);
  }
}

TEST(OpenSource, MissingFileThrows) {
  EXPECT_THROW(io::open_source("/nonexistent/itchbook/file.gz"), std::runtime_error);
}

}  // namespace
}  // namespace itchbook::itch
