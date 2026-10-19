#pragma once
// Helpers shared by stream tests and property tests: random valid messages,
// a handler that records every delivered payload verbatim, and gzip encoding.

#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include <zlib.h>

#include "itchbook/itch/dispatch.hpp"
#include "itchbook/itch/messages.hpp"
#include "support/writer.hpp"

namespace itchbook::test {

inline std::vector<char> known_types() {
  std::vector<char> t;
  for (std::size_t i = 0; i < 256; ++i) {
    if (itch::kMessageSize[i] != 0) {
      t.push_back(static_cast<char>(i));
    }
  }
  return t;
}

// A message of the given type with correct length and random body bytes.
inline std::vector<std::uint8_t> random_message(char type, std::mt19937_64& rng) {
  std::vector<std::uint8_t> m(itch::kMessageSize[static_cast<std::uint8_t>(type)]);
  m[0] = static_cast<std::uint8_t>(type);
  for (std::size_t i = 1; i < m.size(); ++i) {
    m[i] = static_cast<std::uint8_t>(rng());
  }
  return m;
}

struct Recorder : itch::NullHandler {
  using NullHandler::on;
  std::vector<std::vector<std::uint8_t>> payloads;
  template <class M>
  void on(const M& m) {
    payloads.emplace_back(m.data(), m.data() + M::kSize);
  }
};

// gzip-encode `raw` as one gzip member.
inline std::vector<std::uint8_t> gzip(std::span<const std::uint8_t> raw, int level = 6) {
  z_stream zs{};
  if (deflateInit2(&zs, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
    throw std::runtime_error("deflateInit2");
  }
  std::vector<std::uint8_t> out(deflateBound(&zs, static_cast<uLong>(raw.size())) + 32);
  zs.next_in = const_cast<Bytef*>(raw.data());  // NOLINT(cppcoreguidelines-pro-type-const-cast): zlib API
  zs.avail_in = static_cast<uInt>(raw.size());
  zs.next_out = out.data();
  zs.avail_out = static_cast<uInt>(out.size());
  if (deflate(&zs, Z_FINISH) != Z_STREAM_END) {
    deflateEnd(&zs);
    throw std::runtime_error("deflate");
  }
  out.resize(zs.total_out);
  deflateEnd(&zs);
  return out;
}

}  // namespace itchbook::test
