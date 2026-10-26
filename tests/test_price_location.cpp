// price_location: every region of the quote, half-tick mids, and the undefined cases.

#include <gtest/gtest.h>

#include "itchbook/book/price_location.hpp"

namespace itchbook::book {
namespace {

using PL = PriceLocation;

TEST(PriceLocation, EveryRegionOfAOneCentQuote) {
  // bid 10.00, ask 10.02, mid 10.01
  const Price b = 10'0000, a = 10'0200;
  EXPECT_EQ(price_location(9'9900, b, a), PL::BelowBid);
  EXPECT_EQ(price_location(b, b, a), PL::AtBid);
  EXPECT_EQ(price_location(10'0050, b, a), PL::InsideBelowMid);
  EXPECT_EQ(price_location(10'0100, b, a), PL::AtMid);
  EXPECT_EQ(price_location(10'0150, b, a), PL::InsideAboveMid);
  EXPECT_EQ(price_location(a, b, a), PL::AtAsk);
  EXPECT_EQ(price_location(10'0300, b, a), PL::AboveAsk);
}

TEST(PriceLocation, HalfTickMidOnAOneTickSpread) {
  // bid 10.00, ask 10.01: mid 10.005 is exact via 2*px == bid + ask
  const Price b = 10'0000, a = 10'0100;
  EXPECT_EQ(price_location(10'0050, b, a), PL::AtMid);
  EXPECT_EQ(price_location(10'0049, b, a), PL::InsideBelowMid);
  EXPECT_EQ(price_location(10'0051, b, a), PL::InsideAboveMid);
}

TEST(PriceLocation, UndefinedWhenLockedCrossedOrOneSided) {
  EXPECT_EQ(price_location(10'0000, 10'0000, 10'0000), PL::MidUndefined);  // locked
  EXPECT_EQ(price_location(10'0000, 10'0100, 10'0000), PL::MidUndefined);  // crossed
  EXPECT_EQ(price_location(10'0000, 0, 10'0100), PL::MidUndefined);        // no bid
  EXPECT_EQ(price_location(10'0000, 10'0000, 0), PL::MidUndefined);        // no ask
}

TEST(PriceLocation, NoOverflowAtTheTopOfThePriceRange) {
  const Price b = 0xFFFF'FFF0U, a = 0xFFFF'FFFEU;  // sum exceeds 32 bits
  EXPECT_EQ(price_location(0xFFFF'FFF7U, b, a), PL::AtMid);
}

TEST(PriceLocation, Names) {
  EXPECT_EQ(to_string(PL::AtMid), "at_mid");
  EXPECT_EQ(to_string(PL::MidUndefined), "mid_undefined");
}

}  // namespace
}  // namespace itchbook::book
