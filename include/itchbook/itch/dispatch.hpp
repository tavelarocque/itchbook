#pragma once
// Type dispatch: one switch on the type byte, a strict length check, then a
// call to the handler's overload for that message view.
//
// A handler is any object with `on(const M&)` overloads. Deriving from
// NullHandler and writing `using NullHandler::on;` gives a no-op default for
// every type the handler does not care about, and lets the compiler drop
// those branches entirely.

#include <cstddef>
#include <cstdint>

#include "itchbook/itch/messages.hpp"

namespace itchbook::itch {

enum class DecodeStatus : std::uint8_t {
  ok,
  unknown_type,  // type byte is not an ITCH 5.0 message
  bad_length,    // known type, but frame length differs from the spec length
};

struct NullHandler {
  template <class M>
  void on(const M& /*msg*/) noexcept {}
};

namespace detail {
template <class M, class Handler>
inline DecodeStatus deliver(const std::uint8_t* p, std::size_t len, Handler& h) {
  if (len != M::kSize) [[unlikely]] {
    return DecodeStatus::bad_length;
  }
  h.on(M{p});
  return DecodeStatus::ok;
}
}  // namespace detail

// `p` points at the type byte, `len` is the frame length (> 0).
template <class Handler>
inline DecodeStatus decode(const std::uint8_t* p, std::size_t len, Handler& h) {
  using detail::deliver;
  // Ordered roughly by frequency in a full trading day.
  switch (static_cast<char>(p[0])) {
    case 'A': return deliver<AddOrder>(p, len, h);
    case 'D': return deliver<OrderDelete>(p, len, h);
    case 'U': return deliver<OrderReplace>(p, len, h);
    case 'E': return deliver<OrderExecuted>(p, len, h);
    case 'X': return deliver<OrderCancel>(p, len, h);
    case 'I': return deliver<Noii>(p, len, h);
    case 'F': return deliver<AddOrderMpid>(p, len, h);
    case 'P': return deliver<Trade>(p, len, h);
    case 'L': return deliver<MarketParticipantPosition>(p, len, h);
    case 'C': return deliver<OrderExecutedWithPrice>(p, len, h);
    case 'Q': return deliver<CrossTrade>(p, len, h);
    case 'Y': return deliver<RegShoRestriction>(p, len, h);
    case 'H': return deliver<StockTradingAction>(p, len, h);
    case 'R': return deliver<StockDirectory>(p, len, h);
    case 'S': return deliver<SystemEvent>(p, len, h);
    case 'J': return deliver<LuldAuctionCollar>(p, len, h);
    case 'K': return deliver<IpoQuotingPeriod>(p, len, h);
    case 'V': return deliver<MwcbDeclineLevel>(p, len, h);
    case 'W': return deliver<MwcbStatus>(p, len, h);
    case 'h': return deliver<OperationalHalt>(p, len, h);
    case 'B': return deliver<BrokenTrade>(p, len, h);
    case 'N': return deliver<RetailInterest>(p, len, h);
    case 'O': return deliver<DlcrPriceDiscovery>(p, len, h);
    default: return DecodeStatus::unknown_type;
  }
}

}  // namespace itchbook::itch
