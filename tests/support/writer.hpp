#pragma once
// Test-side ITCH encoder. It appends fields one after another in spec order
// (like a struct format string), and never uses the decoder's offsets, so a
// round trip through it checks the decoder's offset table against an
// independent sequential layout.

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace itchbook::test {

class Writer {
 public:
  Writer& u8(std::uint8_t v) {
    b_.push_back(v);
    return *this;
  }
  Writer& ch(char c) { return u8(static_cast<std::uint8_t>(c)); }
  Writer& be(std::uint64_t v, int bytes) {
    for (int i = bytes - 1; i >= 0; --i) {
      b_.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
    }
    return *this;
  }
  Writer& u16(std::uint16_t v) { return be(v, 2); }
  Writer& u32(std::uint32_t v) { return be(v, 4); }
  Writer& u48(std::uint64_t v) { return be(v, 6); }
  Writer& u64(std::uint64_t v) { return be(v, 8); }
  // Left-justified, space-padded ASCII of exactly n bytes.
  Writer& alpha(std::string_view s, std::size_t n) {
    if (s.size() > n) {
      throw std::invalid_argument("alpha field too long");
    }
    for (std::size_t i = 0; i < n; ++i) {
      ch(i < s.size() ? s[i] : ' ');
    }
    return *this;
  }
  // Common 11-byte header.
  Writer& header(char type, std::uint16_t locate, std::uint16_t tracking, std::uint64_t ts) {
    return ch(type).u16(locate).u16(tracking).u48(ts);
  }

  [[nodiscard]] const std::vector<std::uint8_t>& bytes() const { return b_; }
  [[nodiscard]] std::size_t size() const { return b_.size(); }

 private:
  std::vector<std::uint8_t> b_;
};

// Append one BinaryFILE frame (2-byte big-endian length + payload).
inline void append_frame(std::vector<std::uint8_t>& out, std::span<const std::uint8_t> payload) {
  if (payload.size() > 0xFFFF) {
    throw std::invalid_argument("frame too long");
  }
  out.push_back(static_cast<std::uint8_t>(payload.size() >> 8U));
  out.push_back(static_cast<std::uint8_t>(payload.size() & 0xFFU));
  out.insert(out.end(), payload.begin(), payload.end());
}

}  // namespace itchbook::test
