#pragma once
// Frame splitter: turns a byte source into ITCH messages delivered to a
// handler, using one fixed-size buffer.
//
// File framing (NASDAQ BinaryFILE): each message is preceded by a 2-byte
// big-endian length that counts the payload only. A zero-length frame marks
// end of session. A frame can straddle two reads, so after parsing every
// complete frame in the buffer the incomplete tail (at most 2 + 65535 bytes)
// is moved to the front and the next read appends after it. Messages are
// handed to the handler as views into this buffer and are never copied.

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/dispatch.hpp"
#include "itchbook/itch/endian.hpp"

namespace itchbook::itch {

struct StreamLimits {
  // Stop after this many successfully decoded messages (development runs).
  std::uint64_t max_messages = std::numeric_limits<std::uint64_t>::max();
};

struct StreamStats {
  std::uint64_t frames = 0;                 // every length-prefixed frame, including zero-length ones
  std::uint64_t messages = 0;               // frames decoded and delivered to the handler
  std::uint64_t bytes = 0;                  // bytes consumed, length prefixes included
  std::uint64_t end_of_session_frames = 0;  // zero-length frames
  std::uint64_t unknown_type = 0;
  std::uint64_t bad_length = 0;
  std::uint64_t truncated_tail_bytes = 0;  // bytes of an incomplete final frame
  bool source_truncated = false;           // compressed input ended mid-stream
  bool stopped_at_limit = false;

  [[nodiscard]] bool clean() const noexcept {
    return unknown_type == 0 && bad_length == 0 && truncated_tail_bytes == 0 && !source_truncated;
  }
};

inline constexpr std::size_t kMaxFrameBytes = 2 + std::numeric_limits<std::uint16_t>::max();

// Optional software prefetching. A handler that defines
//   void prefetch(const std::uint8_t* payload, std::size_t len, int stage);
// is shown each frame twice before it is decoded: at stage 0 when it is
// kFarAhead frames ahead (start loading what can be located from the message
// alone, e.g. a hash bucket) and at stage 1 when it is kNearAhead frames ahead
// (follow what stage 0 brought into cache, e.g. the order slot the bucket
// names). Prefetching never changes behaviour, only how long loads take.
template <class H>
concept Prefetching = requires(H& h, const std::uint8_t* p, std::size_t n) { h.prefetch(p, n, 0); };
inline constexpr int kFarAhead = 16;
inline constexpr int kNearAhead = 8;

template <class Handler>
StreamStats parse_stream(io::ByteSource& src, Handler& h, const StreamLimits& limits = {},
                         std::size_t read_bytes = std::size_t{1} << 20U) {
  StreamStats st;
  // read_bytes of fresh space is always available after compaction, because
  // the leftover tail is shorter than one maximal frame.
  std::vector<std::uint8_t> buf(read_bytes + kMaxFrameBytes);
  std::size_t begin = 0;
  std::size_t end = 0;

  // Prefetch cursors (used only by Prefetching handlers): byte offsets of the
  // next frame each stage will look at. They run ahead of `begin` within the
  // bytes already in the buffer.
  std::size_t far = 0;   // NOLINT(misc-const-correctness): written only by Prefetching handlers
  std::size_t near = 0;  // NOLINT(misc-const-correctness)
  const auto scout = [&](std::size_t& cur, int stage) {
    if (end - cur >= 2) {
      const std::size_t len = load_be16(buf.data() + cur);
      if (end - cur >= 2 + len) {
        if (len != 0) {
          if constexpr (Prefetching<Handler>) {
            h.prefetch(buf.data() + cur + 2, len, stage);
          }
        }
        cur += 2 + len;
      }
    }
  };

  for (;;) {
    if constexpr (Prefetching<Handler>) {
      far = near = begin;
      for (int i = 0; i < kFarAhead; ++i) {
        scout(far, 0);
      }
      for (int i = 0; i < kNearAhead; ++i) {
        scout(near, 1);
      }
    }
    while (end - begin >= 2) {
      if constexpr (Prefetching<Handler>) {
        scout(far, 0);
        scout(near, 1);
      }
      const std::uint8_t* frame = buf.data() + begin;
      const std::size_t len = load_be16(frame);
      if (end - begin < 2 + len) {
        break;
      }
      if (len != 0 && st.messages == limits.max_messages) {
        st.stopped_at_limit = true;
        return st;
      }
      begin += 2 + len;
      ++st.frames;
      st.bytes += 2 + len;
      if (len == 0) {
        ++st.end_of_session_frames;
        continue;
      }
      switch (decode(frame + 2, len, h)) {
        case DecodeStatus::ok: ++st.messages; break;
        case DecodeStatus::unknown_type: ++st.unknown_type; break;
        case DecodeStatus::bad_length: ++st.bad_length; break;
      }
    }

    if (begin > 0) {
      std::memmove(buf.data(), buf.data() + begin, end - begin);
      end -= begin;
      begin = 0;
    }
    const std::size_t n = src.read(buf.data() + end, buf.size() - end);
    if (n == 0) {
      st.truncated_tail_bytes = end;
      st.source_truncated = src.truncated();
      return st;
    }
    end += n;
  }
}

}  // namespace itchbook::itch
