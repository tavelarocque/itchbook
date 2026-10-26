#pragma once
// TapeBuilder: turns the execution-related ITCH messages of a day into
// TapeRecords, reading the displayed quote BEFORE each message is applied.
//
// Header-only and templated on the book engine, so tapes built over FastBook
// and over ReferenceBook can be compared record for record. Book events are
// forwarded to the book unchanged; the builder never filters anything.
//
// Signing: ITCH has no aggressor flag. For E and C the sign is inferred as the
// opposite of the resting order's side (an execution against a resting buy is
// seller-initiated, sign -1). P carries no usable side (always 'B' since
// 14 Jul 2014), so its sign is 0. Q and B have no aggressor.

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>

#include "itchbook/book/types.hpp"
#include "itchbook/tape/record.hpp"

namespace itchbook::tape {

using book::Locate;
using book::OrderRef;
using book::Price;
using book::Shares;
using book::Side;

template <class Book, class Sink>
class TapeBuilder {
 public:
  TapeBuilder(Book& book, Sink sink) : b_(book), sink_(std::move(sink)) { tstate_.fill('?'); }

  // ---- state ----------------------------------------------------------------
  void system_event(char code) {
    if (code == 'Q') phase_ = kRegular;
    if (code == 'M') phase_ = kPost;
  }
  void trading_state(Locate l, char state) { tstate_[l] = state; }
  // 'h': only halts on the NASDAQ market (market code 'Q') are tracked.
  void operational_halt(Locate l, char market, char action) {
    if (market == 'Q') halted_[l] = action == 'H';
  }

  // ---- book events, forwarded ---------------------------------------------
  void add(Locate l, OrderRef ref, char side, Shares n, Price p) { b_.add(l, ref, side, n, p); }
  void cancel(Locate l, OrderRef ref, Shares n) { b_.cancel(l, ref, n); }
  void remove(Locate l, OrderRef ref) { b_.remove(l, ref); }
  void replace(Locate l, OrderRef old_ref, OrderRef new_ref, Shares n, Price p) {
    b_.replace(l, old_ref, new_ref, n, p);
  }

  // ---- records ----------------------------------------------------------------
  // E (type 'E', printable '-', exec_price ignored) and C (type 'C').
  void execute(std::uint32_t seq, std::uint64_t ts, Locate l, OrderRef ref, Shares n, std::uint64_t match, char type,
               char printable, Price exec_price) {
    TapeRecord r = base(seq, ts, type, match);
    r.ref = ref;
    r.printable = type == 'C' ? printable : '-';
    r.shares = n;
    const auto o = lookup(ref);
    const Locate loc = o ? o->locate : l;
    r.locate = loc;
    quote(r, loc);  // before the execution is applied
    if (o) {
      r.side = book::side_char(o->side);
      r.sign = o->side == Side::buy ? std::int8_t{-1} : std::int8_t{1};
      r.at_head = b_.at_head(ref) ? 1 : 0;
      r.price = type == 'C' ? exec_price : o->price;
      if (n > o->shares) r.flags |= kOverExecute;
    } else {
      r.side = '?';
      r.price = type == 'C' ? exec_price : 0;
      r.flags |= kOrphan;
    }
    state(r, loc);
    b_.execute(l, ref, n);
    sink_(r);
  }

  void hidden(std::uint32_t seq, std::uint64_t ts, Locate l, OrderRef ref, char side, Shares n, Price p,
              std::uint64_t match) {
    TapeRecord r = base(seq, ts, 'P', match);
    r.ref = ref;
    r.side = side;
    r.shares = n;
    r.price = p;
    r.locate = l;
    quote(r, l);
    state(r, l);
    sink_(r);
  }

  void cross(std::uint32_t seq, std::uint64_t ts, Locate l, std::uint64_t n, Price p, std::uint64_t match,
             char cross_type) {
    TapeRecord r = base(seq, ts, 'Q', match);
    r.side = cross_type;
    r.shares = sat32(n, r.flags, kSharesSaturated);
    r.price = p;
    r.locate = l;
    quote(r, l);
    state(r, l);
    sink_(r);
  }

  void broken(std::uint32_t seq, std::uint64_t ts, Locate l, std::uint64_t match) {
    TapeRecord r = base(seq, ts, 'B', match);
    r.locate = l;
    quote(r, l);
    state(r, l);
    sink_(r);
  }

  [[nodiscard]] std::uint8_t phase() const noexcept { return phase_; }

 private:
  struct Rest {
    Locate locate;
    Side side;
    Price price;
    Shares shares;
  };

  // FastBook has an O(1) resting(); ReferenceBook only has order().
  std::optional<Rest> lookup(OrderRef ref) const {
    if constexpr (requires(const Book& bk, OrderRef r) { bk.resting(r); }) {
      const auto o = b_.resting(ref);
      if (!o) return std::nullopt;
      return Rest{o->locate, o->side, o->price, o->shares};
    } else {
      const auto o = b_.order(ref);
      if (!o) return std::nullopt;
      return Rest{o->locate, o->side, o->price, o->shares};
    }
  }

  static TapeRecord base(std::uint32_t seq, std::uint64_t ts, char type, std::uint64_t match) {
    TapeRecord r{};
    r.seq = seq;
    r.ts = ts;
    r.type = type;
    r.match = match;
    r.printable = '-';
    r.side = '-';
    return r;
  }

  static std::uint32_t sat32(std::uint64_t v, std::uint8_t& flags, std::uint8_t bit) {
    if (v > std::numeric_limits<std::uint32_t>::max()) {
      flags |= bit;
      return std::numeric_limits<std::uint32_t>::max();
    }
    return static_cast<std::uint32_t>(v);
  }

  void quote(TapeRecord& r, Locate l) const {
    if (const auto bb = b_.best(l, Side::buy)) {
      r.bid_px = bb->price;
      r.bid_sz = sat32(bb->shares, r.flags, kBidSzSaturated);
    }
    if (const auto ba = b_.best(l, Side::sell)) {
      r.ask_px = ba->price;
      r.ask_sz = sat32(ba->shares, r.flags, kAskSzSaturated);
    }
  }

  void state(TapeRecord& r, Locate l) const {
    r.tstate = tstate_[l];
    r.phase = phase_;
    if (halted_[l]) r.flags |= kOperationalHalt;
  }

  Book& b_;
  Sink sink_;
  std::uint8_t phase_ = kPre;
  std::array<char, 65536> tstate_{};
  std::array<bool, 65536> halted_{};
};

}  // namespace itchbook::tape
