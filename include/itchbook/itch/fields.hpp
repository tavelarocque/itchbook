#pragma once
// Field reflection: visit every field of a message as (name, unsigned value).
//
// Used by the fingerprint in itch-stats (which the independent Python
// reference parser reproduces) and by debugging dumps. Every byte of every
// message except the type byte belongs to exactly one field, so a wrong
// offset or width in any accessor changes some field's sum.
//
// Value encoding: integers as decoded; single characters as their byte value;
// multi-byte alpha fields as the big-endian integer of their bytes (so "AAPL    "
// becomes 0x4141504C20202020).

#include <cstdint>
#include <string_view>

#include "itchbook/itch/messages.hpp"

namespace itchbook::itch {

[[nodiscard]] inline std::uint64_t alpha_as_u64(std::string_view s) noexcept {
  std::uint64_t v = 0;
  for (const char c : s) {
    v = (v << 8U) | static_cast<std::uint8_t>(c);
  }
  return v;
}

namespace detail {
[[nodiscard]] inline std::uint64_t val(char c) noexcept { return static_cast<std::uint8_t>(c); }
[[nodiscard]] inline std::uint64_t val(std::string_view s) noexcept { return alpha_as_u64(s); }
[[nodiscard]] inline std::uint64_t val(std::uint64_t v) noexcept { return v; }
[[nodiscard]] inline std::uint64_t val(std::uint32_t v) noexcept { return v; }
[[nodiscard]] inline std::uint64_t val(std::uint16_t v) noexcept { return v; }

template <class M, class F>
inline void header(const M& m, F& f) {
  f("stock_locate", val(m.stock_locate()));
  f("tracking_number", val(m.tracking_number()));
  f("timestamp", val(m.timestamp()));
}
}  // namespace detail

#define ITCHBOOK_FIELD(name) f(#name, detail::val(m.name()))  // NOLINT(cppcoreguidelines-macro-usage)

template <class F>
inline void visit_fields(const SystemEvent& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(event_code);
}
template <class F>
inline void visit_fields(const StockDirectory& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(market_category);
  ITCHBOOK_FIELD(financial_status);
  ITCHBOOK_FIELD(round_lot_size);
  ITCHBOOK_FIELD(round_lots_only);
  ITCHBOOK_FIELD(issue_classification);
  ITCHBOOK_FIELD(issue_sub_type);
  ITCHBOOK_FIELD(authenticity);
  ITCHBOOK_FIELD(short_sale_threshold);
  ITCHBOOK_FIELD(ipo_flag);
  ITCHBOOK_FIELD(luld_tier);
  ITCHBOOK_FIELD(etp_flag);
  ITCHBOOK_FIELD(etp_leverage_factor);
  ITCHBOOK_FIELD(inverse_indicator);
}
template <class F>
inline void visit_fields(const StockTradingAction& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(trading_state);
  ITCHBOOK_FIELD(reserved);
  ITCHBOOK_FIELD(reason);
}
template <class F>
inline void visit_fields(const RegShoRestriction& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(reg_sho_action);
}
template <class F>
inline void visit_fields(const MarketParticipantPosition& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(mpid);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(primary_market_maker);
  ITCHBOOK_FIELD(market_maker_mode);
  ITCHBOOK_FIELD(market_participant_state);
}
template <class F>
inline void visit_fields(const MwcbDeclineLevel& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(level1);
  ITCHBOOK_FIELD(level2);
  ITCHBOOK_FIELD(level3);
}
template <class F>
inline void visit_fields(const MwcbStatus& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(breached_level);
}
template <class F>
inline void visit_fields(const IpoQuotingPeriod& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(release_time);
  ITCHBOOK_FIELD(release_qualifier);
  ITCHBOOK_FIELD(ipo_price);
}
template <class F>
inline void visit_fields(const LuldAuctionCollar& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(reference_price);
  ITCHBOOK_FIELD(upper_price);
  ITCHBOOK_FIELD(lower_price);
  ITCHBOOK_FIELD(extension);
}
template <class F>
inline void visit_fields(const OperationalHalt& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(market_code);
  ITCHBOOK_FIELD(halt_action);
}
template <class F>
inline void visit_fields(const AddOrder& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(side);
  ITCHBOOK_FIELD(shares);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(price);
}
template <class F>
inline void visit_fields(const AddOrderMpid& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(side);
  ITCHBOOK_FIELD(shares);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(price);
  ITCHBOOK_FIELD(attribution);
}
template <class F>
inline void visit_fields(const OrderExecuted& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(executed_shares);
  ITCHBOOK_FIELD(match_number);
}
template <class F>
inline void visit_fields(const OrderExecutedWithPrice& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(executed_shares);
  ITCHBOOK_FIELD(match_number);
  ITCHBOOK_FIELD(printable);
  ITCHBOOK_FIELD(execution_price);
}
template <class F>
inline void visit_fields(const OrderCancel& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(cancelled_shares);
}
template <class F>
inline void visit_fields(const OrderDelete& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
}
template <class F>
inline void visit_fields(const OrderReplace& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(original_order_ref);
  ITCHBOOK_FIELD(new_order_ref);
  ITCHBOOK_FIELD(shares);
  ITCHBOOK_FIELD(price);
}
template <class F>
inline void visit_fields(const Trade& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(order_ref);
  ITCHBOOK_FIELD(side);
  ITCHBOOK_FIELD(shares);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(price);
  ITCHBOOK_FIELD(match_number);
}
template <class F>
inline void visit_fields(const CrossTrade& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(shares);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(cross_price);
  ITCHBOOK_FIELD(match_number);
  ITCHBOOK_FIELD(cross_type);
}
template <class F>
inline void visit_fields(const BrokenTrade& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(match_number);
}
template <class F>
inline void visit_fields(const Noii& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(paired_shares);
  ITCHBOOK_FIELD(imbalance_shares);
  ITCHBOOK_FIELD(imbalance_direction);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(far_price);
  ITCHBOOK_FIELD(near_price);
  ITCHBOOK_FIELD(current_reference_price);
  ITCHBOOK_FIELD(cross_type);
  ITCHBOOK_FIELD(price_variation_indicator);
}
template <class F>
inline void visit_fields(const RetailInterest& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(interest_flag);
}
template <class F>
inline void visit_fields(const DlcrPriceDiscovery& m, F&& f) {
  detail::header(m, f);
  ITCHBOOK_FIELD(stock);
  ITCHBOOK_FIELD(open_eligibility);
  ITCHBOOK_FIELD(min_allowable_price);
  ITCHBOOK_FIELD(max_allowable_price);
  ITCHBOOK_FIELD(near_execution_price);
  ITCHBOOK_FIELD(near_execution_time);
  ITCHBOOK_FIELD(lower_price_range_collar);
  ITCHBOOK_FIELD(upper_price_range_collar);
}

#undef ITCHBOOK_FIELD

}  // namespace itchbook::itch
