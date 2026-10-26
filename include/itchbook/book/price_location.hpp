#pragma once
// Where a trade price sits relative to a displayed bid and ask.
//
// Used to classify hidden (P) executions against the prevailing NASDAQ displayed
// quote. The mid is defined only when both sides exist and bid < ask; locked,
// crossed and one-sided books give MidUndefined. Prices are ITCH Price(4) units;
// the comparison with the mid is done on 2*price vs bid+ask, so a mid that falls
// on a half tick is handled exactly.

#include <cstdint>
#include <string_view>

#include "itchbook/book/types.hpp"

namespace itchbook::book {

enum class PriceLocation : std::uint8_t {
  MidUndefined,
  BelowBid,
  AtBid,
  InsideBelowMid,
  AtMid,
  InsideAboveMid,
  AtAsk,
  AboveAsk,
};

// bid or ask equal to 0 means that side is empty.
[[nodiscard]] constexpr PriceLocation price_location(Price px, Price bid, Price ask) noexcept {
  if (bid == 0 || ask == 0 || bid >= ask) return PriceLocation::MidUndefined;
  if (px < bid) return PriceLocation::BelowBid;
  if (px == bid) return PriceLocation::AtBid;
  if (px > ask) return PriceLocation::AboveAsk;
  if (px == ask) return PriceLocation::AtAsk;
  const std::uint64_t twice = std::uint64_t{px} * 2U;
  const std::uint64_t sum = std::uint64_t{bid} + ask;
  if (twice == sum) return PriceLocation::AtMid;
  return twice < sum ? PriceLocation::InsideBelowMid : PriceLocation::InsideAboveMid;
}

[[nodiscard]] constexpr std::string_view to_string(PriceLocation l) noexcept {
  switch (l) {
    case PriceLocation::MidUndefined: return "mid_undefined";
    case PriceLocation::BelowBid: return "below_bid";
    case PriceLocation::AtBid: return "at_bid";
    case PriceLocation::InsideBelowMid: return "inside_below_mid";
    case PriceLocation::AtMid: return "at_mid";
    case PriceLocation::InsideAboveMid: return "inside_above_mid";
    case PriceLocation::AtAsk: return "at_ask";
    case PriceLocation::AboveAsk: return "above_ask";
  }
  return "?";
}

}  // namespace itchbook::book
