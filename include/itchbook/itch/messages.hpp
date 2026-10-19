#pragma once
// Zero-copy views over NASDAQ TotalView-ITCH 5.0 messages.
//
// Each view holds one pointer to the first byte of a message payload (the
// type byte; the 2-byte length prefix is not included) and decodes a field
// only when its accessor is called. Constructing a view costs nothing and
// copies nothing. A view is valid only while the buffer it points into is
// alive, which for the stream reader means for the duration of the handler
// callback.
//
// Offsets and lengths follow the ITCH 5.0 specification (revision log to
// 2023-04-28). Every message starts with the same 11-byte header:
//   type@0(1) stock_locate@1(2) tracking_number@3(2) timestamp@5(6)
// Integers are unsigned big-endian. Price(4) fields carry 4 implied decimals
// and are returned as raw integers (so 1234500 means $123.4500). Alpha fields
// are left-justified, space-padded ASCII and are returned as string_views of
// their full width, padding included.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

#include "itchbook/itch/endian.hpp"

namespace itchbook::itch {

template <char Type, std::size_t Size>
class MessageBase {
 public:
  static constexpr char kType = Type;
  static constexpr std::size_t kSize = Size;

  explicit constexpr MessageBase(const std::uint8_t* p) noexcept : p_(p) {}

  [[nodiscard]] const std::uint8_t* data() const noexcept { return p_; }
  [[nodiscard]] char type() const noexcept { return ch(0); }
  [[nodiscard]] std::uint16_t stock_locate() const noexcept { return u16(1); }
  [[nodiscard]] std::uint16_t tracking_number() const noexcept { return u16(3); }
  [[nodiscard]] std::uint64_t timestamp() const noexcept { return load_be48(p_ + 5); }

 protected:
  [[nodiscard]] char ch(std::size_t off) const noexcept { return static_cast<char>(p_[off]); }
  [[nodiscard]] std::uint16_t u16(std::size_t off) const noexcept { return load_be16(p_ + off); }
  [[nodiscard]] std::uint32_t u32(std::size_t off) const noexcept { return load_be32(p_ + off); }
  [[nodiscard]] std::uint64_t u64(std::size_t off) const noexcept { return load_be64(p_ + off); }
  [[nodiscard]] std::string_view alpha(std::size_t off, std::size_t n) const noexcept {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): bytes are ASCII
    return {reinterpret_cast<const char*>(p_ + off), n};
  }

 private:
  const std::uint8_t* p_;
};

// ---- System and reference messages ----------------------------------------

struct SystemEvent : MessageBase<'S', 12> {
  using MessageBase::MessageBase;
  [[nodiscard]] char event_code() const noexcept { return ch(11); }
};

struct StockDirectory : MessageBase<'R', 39> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char market_category() const noexcept { return ch(19); }
  [[nodiscard]] char financial_status() const noexcept { return ch(20); }
  [[nodiscard]] std::uint32_t round_lot_size() const noexcept { return u32(21); }
  [[nodiscard]] char round_lots_only() const noexcept { return ch(25); }
  [[nodiscard]] char issue_classification() const noexcept { return ch(26); }
  [[nodiscard]] std::string_view issue_sub_type() const noexcept { return alpha(27, 2); }
  [[nodiscard]] char authenticity() const noexcept { return ch(29); }
  [[nodiscard]] char short_sale_threshold() const noexcept { return ch(30); }
  [[nodiscard]] char ipo_flag() const noexcept { return ch(31); }
  [[nodiscard]] char luld_tier() const noexcept { return ch(32); }
  [[nodiscard]] char etp_flag() const noexcept { return ch(33); }
  [[nodiscard]] std::uint32_t etp_leverage_factor() const noexcept { return u32(34); }
  [[nodiscard]] char inverse_indicator() const noexcept { return ch(38); }
};

struct StockTradingAction : MessageBase<'H', 25> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char trading_state() const noexcept { return ch(19); }
  [[nodiscard]] char reserved() const noexcept { return ch(20); }
  [[nodiscard]] std::string_view reason() const noexcept { return alpha(21, 4); }
};

struct RegShoRestriction : MessageBase<'Y', 20> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char reg_sho_action() const noexcept { return ch(19); }
};

struct MarketParticipantPosition : MessageBase<'L', 26> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view mpid() const noexcept { return alpha(11, 4); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(15, 8); }
  [[nodiscard]] char primary_market_maker() const noexcept { return ch(23); }
  [[nodiscard]] char market_maker_mode() const noexcept { return ch(24); }
  [[nodiscard]] char market_participant_state() const noexcept { return ch(25); }
};

// Price(8): 8 implied decimals.
struct MwcbDeclineLevel : MessageBase<'V', 35> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t level1() const noexcept { return u64(11); }
  [[nodiscard]] std::uint64_t level2() const noexcept { return u64(19); }
  [[nodiscard]] std::uint64_t level3() const noexcept { return u64(27); }
};

struct MwcbStatus : MessageBase<'W', 12> {
  using MessageBase::MessageBase;
  [[nodiscard]] char breached_level() const noexcept { return ch(11); }
};

struct IpoQuotingPeriod : MessageBase<'K', 28> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  // Seconds since midnight, not nanoseconds.
  [[nodiscard]] std::uint32_t release_time() const noexcept { return u32(19); }
  [[nodiscard]] char release_qualifier() const noexcept { return ch(23); }
  [[nodiscard]] std::uint32_t ipo_price() const noexcept { return u32(24); }
};

struct LuldAuctionCollar : MessageBase<'J', 35> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] std::uint32_t reference_price() const noexcept { return u32(19); }
  [[nodiscard]] std::uint32_t upper_price() const noexcept { return u32(23); }
  [[nodiscard]] std::uint32_t lower_price() const noexcept { return u32(27); }
  [[nodiscard]] std::uint32_t extension() const noexcept { return u32(31); }
};

struct OperationalHalt : MessageBase<'h', 21> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char market_code() const noexcept { return ch(19); }
  [[nodiscard]] char halt_action() const noexcept { return ch(20); }
};

// ---- Order book messages ---------------------------------------------------

struct AddOrder : MessageBase<'A', 36> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] char side() const noexcept { return ch(19); }
  [[nodiscard]] std::uint32_t shares() const noexcept { return u32(20); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(24, 8); }
  [[nodiscard]] std::uint32_t price() const noexcept { return u32(32); }
};

struct AddOrderMpid : MessageBase<'F', 40> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] char side() const noexcept { return ch(19); }
  [[nodiscard]] std::uint32_t shares() const noexcept { return u32(20); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(24, 8); }
  [[nodiscard]] std::uint32_t price() const noexcept { return u32(32); }
  [[nodiscard]] std::string_view attribution() const noexcept { return alpha(36, 4); }
};

struct OrderExecuted : MessageBase<'E', 31> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] std::uint32_t executed_shares() const noexcept { return u32(19); }
  [[nodiscard]] std::uint64_t match_number() const noexcept { return u64(23); }
};

struct OrderExecutedWithPrice : MessageBase<'C', 36> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] std::uint32_t executed_shares() const noexcept { return u32(19); }
  [[nodiscard]] std::uint64_t match_number() const noexcept { return u64(23); }
  [[nodiscard]] char printable() const noexcept { return ch(31); }
  [[nodiscard]] std::uint32_t execution_price() const noexcept { return u32(32); }
};

struct OrderCancel : MessageBase<'X', 23> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] std::uint32_t cancelled_shares() const noexcept { return u32(19); }
};

struct OrderDelete : MessageBase<'D', 19> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
};

struct OrderReplace : MessageBase<'U', 35> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t original_order_ref() const noexcept { return u64(11); }
  [[nodiscard]] std::uint64_t new_order_ref() const noexcept { return u64(19); }
  [[nodiscard]] std::uint32_t shares() const noexcept { return u32(27); }
  [[nodiscard]] std::uint32_t price() const noexcept { return u32(31); }
};

// ---- Trade messages --------------------------------------------------------

struct Trade : MessageBase<'P', 44> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t order_ref() const noexcept { return u64(11); }
  [[nodiscard]] char side() const noexcept { return ch(19); }
  [[nodiscard]] std::uint32_t shares() const noexcept { return u32(20); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(24, 8); }
  [[nodiscard]] std::uint32_t price() const noexcept { return u32(32); }
  [[nodiscard]] std::uint64_t match_number() const noexcept { return u64(36); }
};

struct CrossTrade : MessageBase<'Q', 40> {
  using MessageBase::MessageBase;
  // 8 bytes here, unlike every other shares field.
  [[nodiscard]] std::uint64_t shares() const noexcept { return u64(11); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(19, 8); }
  [[nodiscard]] std::uint32_t cross_price() const noexcept { return u32(27); }
  [[nodiscard]] std::uint64_t match_number() const noexcept { return u64(31); }
  [[nodiscard]] char cross_type() const noexcept { return ch(39); }
};

struct BrokenTrade : MessageBase<'B', 19> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t match_number() const noexcept { return u64(11); }
};

// ---- Auction and retail messages --------------------------------------------

struct Noii : MessageBase<'I', 50> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::uint64_t paired_shares() const noexcept { return u64(11); }
  [[nodiscard]] std::uint64_t imbalance_shares() const noexcept { return u64(19); }
  [[nodiscard]] char imbalance_direction() const noexcept { return ch(27); }
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(28, 8); }
  [[nodiscard]] std::uint32_t far_price() const noexcept { return u32(36); }
  [[nodiscard]] std::uint32_t near_price() const noexcept { return u32(40); }
  [[nodiscard]] std::uint32_t current_reference_price() const noexcept { return u32(44); }
  [[nodiscard]] char cross_type() const noexcept { return ch(48); }
  [[nodiscard]] char price_variation_indicator() const noexcept { return ch(49); }
};

struct RetailInterest : MessageBase<'N', 20> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char interest_flag() const noexcept { return ch(19); }
};

struct DlcrPriceDiscovery : MessageBase<'O', 48> {
  using MessageBase::MessageBase;
  [[nodiscard]] std::string_view stock() const noexcept { return alpha(11, 8); }
  [[nodiscard]] char open_eligibility() const noexcept { return ch(19); }
  [[nodiscard]] std::uint32_t min_allowable_price() const noexcept { return u32(20); }
  [[nodiscard]] std::uint32_t max_allowable_price() const noexcept { return u32(24); }
  [[nodiscard]] std::uint32_t near_execution_price() const noexcept { return u32(28); }
  [[nodiscard]] std::uint64_t near_execution_time() const noexcept { return u64(32); }
  [[nodiscard]] std::uint32_t lower_price_range_collar() const noexcept { return u32(40); }
  [[nodiscard]] std::uint32_t upper_price_range_collar() const noexcept { return u32(44); }
};

// ---- Type list -------------------------------------------------------------

template <class... Ms>
struct TypeList {};

using AllMessages =
    TypeList<SystemEvent, StockDirectory, StockTradingAction, RegShoRestriction, MarketParticipantPosition,
             MwcbDeclineLevel, MwcbStatus, IpoQuotingPeriod, LuldAuctionCollar, OperationalHalt, AddOrder, AddOrderMpid,
             OrderExecuted, OrderExecutedWithPrice, OrderCancel, OrderDelete, OrderReplace, Trade, CrossTrade,
             BrokenTrade, Noii, RetailInterest, DlcrPriceDiscovery>;

inline constexpr std::size_t kMessageTypeCount = 23;

// Expected total length (type byte included) for each type byte; 0 = unknown.
inline constexpr std::array<std::uint8_t, 256> kMessageSize = [] {
  std::array<std::uint8_t, 256> t{};
  auto set = []<class... Ms>(std::array<std::uint8_t, 256>& a, TypeList<Ms...>) {
    ((a[static_cast<std::uint8_t>(Ms::kType)] = static_cast<std::uint8_t>(Ms::kSize)), ...);
  };
  set(t, AllMessages{});
  return t;
}();

}  // namespace itchbook::itch
