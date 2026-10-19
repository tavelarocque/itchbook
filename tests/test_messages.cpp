// Every message type: encode with the sequential test writer, decode with the
// view, check every accessor. Values use the high bit of every field width so
// sign-extension and truncation bugs show up.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include "itchbook/itch/dispatch.hpp"
#include "itchbook/itch/endian.hpp"
#include "itchbook/itch/messages.hpp"
#include "support/writer.hpp"

namespace itchbook::itch {
namespace {

using test::Writer;

constexpr std::uint16_t kLocate = 0xBEEF;
constexpr std::uint16_t kTracking = 0x8123;
constexpr std::uint64_t kTs = 0xABCDEF012345ULL;  // full 48 bits
constexpr std::uint32_t kU32 = 0xF1E2D3C4U;
constexpr std::uint64_t kU64 = 0xFEDCBA9876543210ULL;

Writer hdr(char type) {
  Writer w;
  w.header(type, kLocate, kTracking, kTs);
  return w;
}

template <class M>
M view(const Writer& w) {
  EXPECT_EQ(w.size(), M::kSize) << "encoded length for type " << M::kType;
  return M{w.bytes().data()};
}

template <class M>
void expect_header(const M& m) {
  EXPECT_EQ(m.type(), M::kType);
  EXPECT_EQ(m.stock_locate(), kLocate);
  EXPECT_EQ(m.tracking_number(), kTracking);
  EXPECT_EQ(m.timestamp(), kTs);
}

TEST(Endian, LoadsAreBigEndianAndUnaligned) {
  const std::vector<std::uint8_t> b = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09};
  EXPECT_EQ(load_be16(b.data() + 1), 0x0102);
  EXPECT_EQ(load_be32(b.data() + 1), 0x01020304U);
  EXPECT_EQ(load_be48(b.data() + 1), 0x010203040506ULL);
  EXPECT_EQ(load_be64(b.data() + 1), 0x0102030405060708ULL);
}

TEST(Messages, SpecLengths) {
  // Total lengths from the ITCH 5.0 specification, type byte included.
  const std::vector<std::pair<char, int>> spec = {{'S', 12}, {'R', 39}, {'H', 25}, {'Y', 20}, {'L', 26}, {'V', 35},
                                                  {'W', 12}, {'K', 28}, {'J', 35}, {'h', 21}, {'A', 36}, {'F', 40},
                                                  {'E', 31}, {'C', 36}, {'X', 23}, {'D', 19}, {'U', 35}, {'P', 44},
                                                  {'Q', 40}, {'B', 19}, {'I', 50}, {'N', 20}, {'O', 48}};
  ASSERT_EQ(spec.size(), kMessageTypeCount);
  int known = 0;
  for (std::size_t t = 0; t < kMessageSize.size(); ++t) {
    known += kMessageSize[t] != 0 ? 1 : 0;
  }
  EXPECT_EQ(known, static_cast<int>(kMessageTypeCount));
  for (const auto& [type, len] : spec) {
    EXPECT_EQ(kMessageSize[static_cast<std::uint8_t>(type)], len) << type;
  }
}

TEST(Messages, SystemEvent) {
  auto w = hdr('S');
  w.ch('O');
  const auto m = view<SystemEvent>(w);
  expect_header(m);
  EXPECT_EQ(m.event_code(), 'O');
}

TEST(Messages, StockDirectory) {
  auto w = hdr('R');
  w.alpha("AAPL", 8).ch('Q').ch('N').u32(kU32).ch('N').ch('C').alpha("Z", 2).ch('P').ch('N').ch(' ').ch('1');
  w.ch('Y').u32(0x80000003U).ch('N');
  const auto m = view<StockDirectory>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "AAPL    ");
  EXPECT_EQ(m.market_category(), 'Q');
  EXPECT_EQ(m.financial_status(), 'N');
  EXPECT_EQ(m.round_lot_size(), kU32);
  EXPECT_EQ(m.round_lots_only(), 'N');
  EXPECT_EQ(m.issue_classification(), 'C');
  EXPECT_EQ(m.issue_sub_type(), "Z ");
  EXPECT_EQ(m.authenticity(), 'P');
  EXPECT_EQ(m.short_sale_threshold(), 'N');
  EXPECT_EQ(m.ipo_flag(), ' ');
  EXPECT_EQ(m.luld_tier(), '1');
  EXPECT_EQ(m.etp_flag(), 'Y');
  EXPECT_EQ(m.etp_leverage_factor(), 0x80000003U);
  EXPECT_EQ(m.inverse_indicator(), 'N');
}

TEST(Messages, StockTradingAction) {
  auto w = hdr('H');
  w.alpha("MSFT", 8).ch('H').ch(' ').alpha("LUDP", 4);
  const auto m = view<StockTradingAction>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "MSFT    ");
  EXPECT_EQ(m.trading_state(), 'H');
  EXPECT_EQ(m.reserved(), ' ');
  EXPECT_EQ(m.reason(), "LUDP");
}

TEST(Messages, RegShoRestriction) {
  auto w = hdr('Y');
  w.alpha("ZVZZT", 8).ch('1');
  const auto m = view<RegShoRestriction>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "ZVZZT   ");
  EXPECT_EQ(m.reg_sho_action(), '1');
}

TEST(Messages, MarketParticipantPosition) {
  auto w = hdr('L');
  w.alpha("GSCO", 4).alpha("SPY", 8).ch('Y').ch('N').ch('A');
  const auto m = view<MarketParticipantPosition>(w);
  expect_header(m);
  EXPECT_EQ(m.mpid(), "GSCO");
  EXPECT_EQ(m.stock(), "SPY     ");
  EXPECT_EQ(m.primary_market_maker(), 'Y');
  EXPECT_EQ(m.market_maker_mode(), 'N');
  EXPECT_EQ(m.market_participant_state(), 'A');
}

TEST(Messages, MwcbDeclineLevel) {
  auto w = hdr('V');
  w.u64(kU64).u64(kU64 - 1).u64(0x8000000000000001ULL);
  const auto m = view<MwcbDeclineLevel>(w);
  expect_header(m);
  EXPECT_EQ(m.level1(), kU64);
  EXPECT_EQ(m.level2(), kU64 - 1);
  EXPECT_EQ(m.level3(), 0x8000000000000001ULL);
}

TEST(Messages, MwcbStatus) {
  auto w = hdr('W');
  w.ch('2');
  const auto m = view<MwcbStatus>(w);
  expect_header(m);
  EXPECT_EQ(m.breached_level(), '2');
}

TEST(Messages, IpoQuotingPeriod) {
  auto w = hdr('K');
  w.alpha("NEWCO", 8).u32(34200).ch('A').u32(kU32);
  const auto m = view<IpoQuotingPeriod>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "NEWCO   ");
  EXPECT_EQ(m.release_time(), 34200U);
  EXPECT_EQ(m.release_qualifier(), 'A');
  EXPECT_EQ(m.ipo_price(), kU32);
}

TEST(Messages, LuldAuctionCollar) {
  auto w = hdr('J');
  w.alpha("TSLA", 8).u32(1).u32(0x80000000U).u32(kU32).u32(7);
  const auto m = view<LuldAuctionCollar>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "TSLA    ");
  EXPECT_EQ(m.reference_price(), 1U);
  EXPECT_EQ(m.upper_price(), 0x80000000U);
  EXPECT_EQ(m.lower_price(), kU32);
  EXPECT_EQ(m.extension(), 7U);
}

TEST(Messages, OperationalHalt) {
  auto w = hdr('h');
  w.alpha("QQQ", 8).ch('Q').ch('H');
  const auto m = view<OperationalHalt>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "QQQ     ");
  EXPECT_EQ(m.market_code(), 'Q');
  EXPECT_EQ(m.halt_action(), 'H');
}

TEST(Messages, AddOrder) {
  auto w = hdr('A');
  w.u64(kU64).ch('B').u32(kU32).alpha("AAPL", 8).u32(1234500);
  const auto m = view<AddOrder>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
  EXPECT_EQ(m.side(), 'B');
  EXPECT_EQ(m.shares(), kU32);
  EXPECT_EQ(m.stock(), "AAPL    ");
  EXPECT_EQ(m.price(), 1234500U);
}

TEST(Messages, AddOrderMpid) {
  auto w = hdr('F');
  w.u64(kU64).ch('S').u32(100).alpha("AMZN", 8).u32(0x77359400U).alpha("NSDQ", 4);
  const auto m = view<AddOrderMpid>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
  EXPECT_EQ(m.side(), 'S');
  EXPECT_EQ(m.shares(), 100U);
  EXPECT_EQ(m.stock(), "AMZN    ");
  EXPECT_EQ(m.price(), 0x77359400U);  // spec maximum, 200000.0000
  EXPECT_EQ(m.attribution(), "NSDQ");
}

TEST(Messages, OrderExecuted) {
  auto w = hdr('E');
  w.u64(kU64).u32(kU32).u64(kU64 ^ 0xFFU);
  const auto m = view<OrderExecuted>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
  EXPECT_EQ(m.executed_shares(), kU32);
  EXPECT_EQ(m.match_number(), kU64 ^ 0xFFU);
}

TEST(Messages, OrderExecutedWithPrice) {
  auto w = hdr('C');
  w.u64(kU64).u32(kU32).u64(42).ch('N').u32(0x80000001U);
  const auto m = view<OrderExecutedWithPrice>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
  EXPECT_EQ(m.executed_shares(), kU32);
  EXPECT_EQ(m.match_number(), 42U);
  EXPECT_EQ(m.printable(), 'N');
  EXPECT_EQ(m.execution_price(), 0x80000001U);
}

TEST(Messages, OrderCancel) {
  auto w = hdr('X');
  w.u64(kU64).u32(kU32);
  const auto m = view<OrderCancel>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
  EXPECT_EQ(m.cancelled_shares(), kU32);
}

TEST(Messages, OrderDelete) {
  auto w = hdr('D');
  w.u64(kU64);
  const auto m = view<OrderDelete>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), kU64);
}

TEST(Messages, OrderReplace) {
  auto w = hdr('U');
  w.u64(kU64).u64(kU64 - 7).u32(kU32).u32(0x80000000U);
  const auto m = view<OrderReplace>(w);
  expect_header(m);
  EXPECT_EQ(m.original_order_ref(), kU64);
  EXPECT_EQ(m.new_order_ref(), kU64 - 7);
  EXPECT_EQ(m.shares(), kU32);
  EXPECT_EQ(m.price(), 0x80000000U);
}

TEST(Messages, Trade) {
  auto w = hdr('P');
  w.u64(0).ch('B').u32(kU32).alpha("GOOG", 8).u32(99).u64(kU64);
  const auto m = view<Trade>(w);
  expect_header(m);
  EXPECT_EQ(m.order_ref(), 0U);
  EXPECT_EQ(m.side(), 'B');
  EXPECT_EQ(m.shares(), kU32);
  EXPECT_EQ(m.stock(), "GOOG    ");
  EXPECT_EQ(m.price(), 99U);
  EXPECT_EQ(m.match_number(), kU64);
}

TEST(Messages, CrossTrade) {
  auto w = hdr('Q');
  w.u64(kU64).alpha("IBM", 8).u32(kU32).u64(5).ch('C');
  const auto m = view<CrossTrade>(w);
  expect_header(m);
  EXPECT_EQ(m.shares(), kU64);  // 8-byte shares field
  EXPECT_EQ(m.stock(), "IBM     ");
  EXPECT_EQ(m.cross_price(), kU32);
  EXPECT_EQ(m.match_number(), 5U);
  EXPECT_EQ(m.cross_type(), 'C');
}

TEST(Messages, BrokenTrade) {
  auto w = hdr('B');
  w.u64(kU64);
  const auto m = view<BrokenTrade>(w);
  expect_header(m);
  EXPECT_EQ(m.match_number(), kU64);
}

TEST(Messages, Noii) {
  auto w = hdr('I');
  w.u64(kU64).u64(kU64 >> 1U).ch('S').alpha("NVDA", 8).u32(1).u32(2).u32(kU32).ch('O').ch('A');
  const auto m = view<Noii>(w);
  expect_header(m);
  EXPECT_EQ(m.paired_shares(), kU64);
  EXPECT_EQ(m.imbalance_shares(), kU64 >> 1U);
  EXPECT_EQ(m.imbalance_direction(), 'S');
  EXPECT_EQ(m.stock(), "NVDA    ");
  EXPECT_EQ(m.far_price(), 1U);
  EXPECT_EQ(m.near_price(), 2U);
  EXPECT_EQ(m.current_reference_price(), kU32);
  EXPECT_EQ(m.cross_type(), 'O');
  EXPECT_EQ(m.price_variation_indicator(), 'A');
}

TEST(Messages, RetailInterest) {
  auto w = hdr('N');
  w.alpha("RPI", 8).ch('A');
  const auto m = view<RetailInterest>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "RPI     ");
  EXPECT_EQ(m.interest_flag(), 'A');
}

TEST(Messages, DlcrPriceDiscovery) {
  auto w = hdr('O');
  w.alpha("DLCR", 8).ch('Y').u32(1).u32(2).u32(3).u64(kU64).u32(0x80000000U).u32(kU32);
  const auto m = view<DlcrPriceDiscovery>(w);
  expect_header(m);
  EXPECT_EQ(m.stock(), "DLCR    ");
  EXPECT_EQ(m.open_eligibility(), 'Y');
  EXPECT_EQ(m.min_allowable_price(), 1U);
  EXPECT_EQ(m.max_allowable_price(), 2U);
  EXPECT_EQ(m.near_execution_price(), 3U);
  EXPECT_EQ(m.near_execution_time(), kU64);
  EXPECT_EQ(m.lower_price_range_collar(), 0x80000000U);
  EXPECT_EQ(m.upper_price_range_collar(), kU32);
}

// ---- Dispatch --------------------------------------------------------------

struct TypeRecorder : NullHandler {
  using NullHandler::on;
  std::string seen;
  template <class M>
  void on(const M& m) {
    seen.push_back(m.type());
  }
};

TEST(Dispatch, EveryKnownTypeReachesItsOverload) {
  for (std::size_t t = 0; t < 256; ++t) {
    std::vector<std::uint8_t> payload(kMessageSize[t] == 0 ? 20 : kMessageSize[t], 0);
    payload[0] = static_cast<std::uint8_t>(t);
    TypeRecorder rec;
    const auto status = decode(payload.data(), payload.size(), rec);
    if (kMessageSize[t] == 0) {
      EXPECT_EQ(status, DecodeStatus::unknown_type) << t;
      EXPECT_TRUE(rec.seen.empty());
    } else {
      EXPECT_EQ(status, DecodeStatus::ok) << t;
      EXPECT_EQ(rec.seen, std::string(1, static_cast<char>(t)));
    }
  }
}

TEST(Dispatch, WrongLengthIsRejectedNotDelivered) {
  for (std::size_t t = 0; t < 256; ++t) {
    if (kMessageSize[t] == 0) {
      continue;
    }
    for (const std::size_t len : {std::size_t{kMessageSize[t]} - 1, std::size_t{kMessageSize[t]} + 1}) {
      std::vector<std::uint8_t> payload(len, 0);
      payload[0] = static_cast<std::uint8_t>(t);
      TypeRecorder rec;
      EXPECT_EQ(decode(payload.data(), payload.size(), rec), DecodeStatus::bad_length) << t << " len " << len;
      EXPECT_TRUE(rec.seen.empty());
    }
  }
}

}  // namespace
}  // namespace itchbook::itch
