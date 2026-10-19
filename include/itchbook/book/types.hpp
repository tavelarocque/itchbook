#pragma once

#include <cstdint>
#include <initializer_list>

namespace itchbook::book {

using OrderRef = std::uint64_t;
using Locate = std::uint16_t;
using Price = std::uint32_t;  // ITCH Price(4): 4 implied decimals
using Shares = std::uint32_t;

enum class Side : std::uint8_t { buy = 0, sell = 1 };

[[nodiscard]] constexpr char side_char(Side s) noexcept { return s == Side::buy ? 'B' : 'S'; }

struct LevelInfo {
  Price price = 0;
  std::uint64_t shares = 0;  // sum over the level's orders (can exceed 32 bits)
  std::uint32_t orders = 0;
  friend bool operator==(const LevelInfo&, const LevelInfo&) = default;
};

// SplitMix64 finaliser: a cheap, well-mixed 64-bit hash. The independent
// Python reference implements the same function, so book checksums can be
// compared across implementations.
[[nodiscard]] constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30U)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27U)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31U);
}

[[nodiscard]] constexpr std::uint64_t mix(std::initializer_list<std::uint64_t> vs) noexcept {
  std::uint64_t h = 0;
  for (const std::uint64_t v : vs) {
    h = splitmix64(h ^ v);
  }
  return h;
}

}  // namespace itchbook::book
