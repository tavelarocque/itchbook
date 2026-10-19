#pragma once
// Unaligned big-endian loads. ITCH frames sit at arbitrary byte offsets in the
// input buffer, so every load goes through memcpy (which the compiler turns
// into a single unaligned load) followed by a byte swap on little-endian hosts.

#include <bit>
#include <cstdint>
#include <cstring>

namespace itchbook::itch {

static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big,
              "mixed-endian hosts are not supported");

[[nodiscard]] inline std::uint16_t load_be16(const std::uint8_t* p) noexcept {
  std::uint16_t v;
  std::memcpy(&v, p, sizeof v);
  if constexpr (std::endian::native == std::endian::little) {
    v = __builtin_bswap16(v);
  }
  return v;
}

[[nodiscard]] inline std::uint32_t load_be32(const std::uint8_t* p) noexcept {
  std::uint32_t v;
  std::memcpy(&v, p, sizeof v);
  if constexpr (std::endian::native == std::endian::little) {
    v = __builtin_bswap32(v);
  }
  return v;
}

[[nodiscard]] inline std::uint64_t load_be64(const std::uint8_t* p) noexcept {
  std::uint64_t v;
  std::memcpy(&v, p, sizeof v);
  if constexpr (std::endian::native == std::endian::little) {
    v = __builtin_bswap64(v);
  }
  return v;
}

// 6-byte ITCH timestamp: nanoseconds since midnight, fits in 48 bits.
[[nodiscard]] inline std::uint64_t load_be48(const std::uint8_t* p) noexcept {
  return (static_cast<std::uint64_t>(load_be16(p)) << 32) | load_be32(p + 2);
}

}  // namespace itchbook::itch
