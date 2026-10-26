// TapeBuilder: unit tests of every field rule, then two properties on random
// event streams: (1) the FastBook tape equals a record-by-record expectation
// computed from a brute-force flat order list written here, independently of
// both books; (2) FastBook and ReferenceBook tapes are byte-identical.

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/reference_book.hpp"
#include "itchbook/tape/builder.hpp"

namespace itchbook::tape {
namespace {

using book::FastBook;
using book::ReferenceBook;

struct Collect {
  std::vector<TapeRecord>* out;
  void operator()(const TapeRecord& r) const { out->push_back(r); }
};

template <class Book>
struct Fx {
  Book book;
  std::vector<TapeRecord> recs;
  TapeBuilder<Book, Collect> tb{book, Collect{&recs}};
};

constexpr Price kP = 10'0000;  // $10.00

TEST(Tape, SignIsTheOppositeOfTheRestingSide) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.add(1, 11, 'S', 100, kP + 100);
  f.tb.execute(5, 1000, 1, 10, 40, 77, 'E', '-', 0);  // against a resting buy: seller-initiated
  f.tb.execute(6, 1000, 1, 11, 40, 78, 'E', '-', 0);  // against a resting sell: buyer-initiated
  ASSERT_EQ(f.recs.size(), 2U);
  EXPECT_EQ(f.recs[0].side, 'B');
  EXPECT_EQ(f.recs[0].sign, -1);
  EXPECT_EQ(f.recs[1].side, 'S');
  EXPECT_EQ(f.recs[1].sign, 1);
  EXPECT_EQ(f.recs[0].seq, 5U);
  EXPECT_EQ(f.recs[0].ts, 1000U);
  EXPECT_EQ(f.recs[0].match, 77U);
  EXPECT_EQ(f.recs[0].ref, 10U);
  EXPECT_EQ(f.recs[0].shares, 40U);
  EXPECT_EQ(f.recs[0].price, kP);  // E: the resting order's price
  EXPECT_EQ(f.recs[0].printable, '-');
  EXPECT_EQ(f.recs[0].type, 'E');
}

TEST(Tape, QuoteIsReadBeforeTheExecution) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'S', 100, kP + 100);
  f.tb.add(1, 11, 'S', 300, kP + 200);
  f.tb.add(1, 12, 'B', 50, kP);
  f.tb.execute(1, 1, 1, 10, 100, 1, 'E', '-', 0);  // empties the best ask level
  f.tb.execute(2, 1, 1, 11, 10, 2, 'E', '-', 0);
  ASSERT_EQ(f.recs.size(), 2U);
  EXPECT_EQ(f.recs[0].ask_px, kP + 100);
  EXPECT_EQ(f.recs[0].ask_sz, 100U);
  EXPECT_EQ(f.recs[0].bid_px, kP);
  EXPECT_EQ(f.recs[0].bid_sz, 50U);
  EXPECT_EQ(f.recs[1].ask_px, kP + 200);
  EXPECT_EQ(f.recs[1].ask_sz, 300U);
}

TEST(Tape, EmptySideIsZero) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.execute(1, 1, 1, 10, 100, 1, 'E', '-', 0);
  EXPECT_EQ(f.recs[0].ask_px, 0U);
  EXPECT_EQ(f.recs[0].ask_sz, 0U);
  EXPECT_EQ(f.recs[0].bid_px, kP);
}

TEST(Tape, CCarriesExecutionPriceAndPrintableFlag) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.execute(1, 1, 1, 10, 30, 9, 'C', 'N', kP + 500);
  f.tb.execute(2, 1, 1, 10, 30, 10, 'C', 'Y', kP + 300);
  EXPECT_EQ(f.recs[0].price, kP + 500);
  EXPECT_EQ(f.recs[0].printable, 'N');
  EXPECT_EQ(f.recs[1].printable, 'Y');
  EXPECT_EQ(f.recs[0].sign, -1);  // still signed; filtering non-printables is downstream
}

TEST(Tape, OrphanAndOverExecuteAreFlagged) {
  Fx<FastBook> f;
  f.tb.execute(1, 1, 3, 999, 10, 1, 'E', '-', 0);
  EXPECT_EQ(f.recs[0].side, '?');
  EXPECT_EQ(f.recs[0].sign, 0);
  EXPECT_EQ(f.recs[0].locate, 3);
  EXPECT_EQ(f.recs[0].price, 0U);
  EXPECT_NE(f.recs[0].flags & kOrphan, 0);
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.execute(2, 1, 1, 10, 150, 2, 'E', '-', 0);
  EXPECT_EQ(f.recs[1].shares, 150U);  // as sent, not clamped
  EXPECT_NE(f.recs[1].flags & kOverExecute, 0);
  EXPECT_EQ(f.recs[1].flags & kOrphan, 0);
}

TEST(Tape, AtHeadIsTheQueuePositionBeforeTheFill) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.add(1, 11, 'B', 100, kP);
  f.tb.execute(1, 1, 1, 11, 10, 1, 'E', '-', 0);   // behind 10
  f.tb.execute(2, 1, 1, 10, 100, 2, 'E', '-', 0);  // at head, removes it
  f.tb.execute(3, 1, 1, 11, 10, 3, 'E', '-', 0);   // now at head
  EXPECT_EQ(f.recs[0].at_head, 0);
  EXPECT_EQ(f.recs[1].at_head, 1);
  EXPECT_EQ(f.recs[2].at_head, 1);
}

TEST(Tape, LocateAndQuoteComeFromTheRestingOrder) {
  Fx<FastBook> f;
  f.tb.add(2, 10, 'B', 100, kP);
  f.tb.add(1, 20, 'B', 7, kP - 100);
  f.tb.execute(1, 1, 1, 10, 10, 1, 'E', '-', 0);  // header locate 1, order on locate 2
  EXPECT_EQ(f.recs[0].locate, 2);
  EXPECT_EQ(f.recs[0].bid_px, kP);
  EXPECT_EQ(f.recs[0].bid_sz, 100U);
}

TEST(Tape, PhaseTradingStateAndOperationalHalt) {
  Fx<FastBook> f;
  f.tb.hidden(1, 1, 1, 0, 'B', 5, kP, 1);
  f.tb.trading_state(1, 'H');
  f.tb.system_event('Q');
  f.tb.hidden(2, 2, 1, 0, 'B', 5, kP, 2);
  f.tb.operational_halt(1, 'B', 'H');  // BX halt: not tracked
  f.tb.hidden(3, 3, 1, 0, 'B', 5, kP, 3);
  f.tb.operational_halt(1, 'Q', 'H');
  f.tb.trading_state(1, 'T');
  f.tb.hidden(4, 4, 1, 0, 'B', 5, kP, 4);
  f.tb.operational_halt(1, 'Q', 'T');
  f.tb.system_event('M');
  f.tb.hidden(5, 5, 1, 0, 'B', 5, kP, 5);
  EXPECT_EQ(f.recs[0].phase, kPre);
  EXPECT_EQ(f.recs[0].tstate, '?');
  EXPECT_EQ(f.recs[1].phase, kRegular);
  EXPECT_EQ(f.recs[1].tstate, 'H');
  EXPECT_EQ(f.recs[2].flags & kOperationalHalt, 0);
  EXPECT_NE(f.recs[3].flags & kOperationalHalt, 0);
  EXPECT_EQ(f.recs[3].tstate, 'T');
  EXPECT_EQ(f.recs[4].phase, kPost);
  EXPECT_EQ(f.recs[4].flags & kOperationalHalt, 0);
}

TEST(Tape, HiddenCrossAndBrokenRecords) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 100, kP);
  f.tb.hidden(1, 1, 1, 0, 'B', 25, kP + 50, 11);
  f.tb.cross(2, 2, 1, std::uint64_t{5'000'000'000}, kP + 10, 12, 'O');
  f.tb.cross(3, 3, 1, 1000, kP + 10, 13, 'C');
  f.tb.broken(4, 4, 1, 11);
  EXPECT_EQ(f.recs[0].type, 'P');
  EXPECT_EQ(f.recs[0].side, 'B');
  EXPECT_EQ(f.recs[0].sign, 0);
  EXPECT_EQ(f.recs[0].price, kP + 50);
  EXPECT_EQ(f.recs[0].bid_px, kP);
  EXPECT_EQ(f.recs[1].type, 'Q');
  EXPECT_EQ(f.recs[1].side, 'O');
  EXPECT_EQ(f.recs[1].shares, std::numeric_limits<std::uint32_t>::max());
  EXPECT_NE(f.recs[1].flags & kSharesSaturated, 0);
  EXPECT_EQ(f.recs[2].shares, 1000U);
  EXPECT_EQ(f.recs[2].flags & kSharesSaturated, 0);
  EXPECT_EQ(f.recs[3].type, 'B');
  EXPECT_EQ(f.recs[3].match, 11U);
  EXPECT_EQ(f.recs[3].price, 0U);
  EXPECT_EQ(f.recs[3].side, '-');
  EXPECT_EQ(f.book.live_orders(), 1U);  // none of P, Q, B touched the book
}

TEST(Tape, LevelSizeSaturates) {
  Fx<FastBook> f;
  f.tb.add(1, 10, 'B', 3'000'000'000U, kP);
  f.tb.add(1, 11, 'B', 3'000'000'000U, kP);
  f.tb.hidden(1, 1, 1, 0, 'B', 1, kP, 1);
  EXPECT_EQ(f.recs[0].bid_sz, std::numeric_limits<std::uint32_t>::max());
  EXPECT_NE(f.recs[0].flags & kBidSzSaturated, 0);
}

// ---- brute-force model -----------------------------------------------------------

struct Model {
  struct O {
    Locate l;
    Side s;
    Price p;
    Shares n;
    std::uint64_t arrival;
  };
  std::map<OrderRef, O> live;
  std::uint64_t clock = 0;
  std::uint8_t phase = kPre;
  std::array<char, 8> tstate{'?', '?', '?', '?', '?', '?', '?', '?'};
  std::array<bool, 8> halted{};

  void add(Locate l, OrderRef r, char sd, Shares n, Price p) {
    live[r] = {l, sd == 'B' ? Side::buy : Side::sell, p, n, ++clock};
  }
  void reduce(OrderRef r, Shares n) {
    auto it = live.find(r);
    if (it == live.end()) return;
    if (n >= it->second.n) {
      live.erase(it);
    } else {
      it->second.n -= n;
    }
  }
  void remove(OrderRef r) { live.erase(r); }
  void replace(OrderRef old_ref, OrderRef new_ref, Shares n, Price p) {
    auto it = live.find(old_ref);
    if (it == live.end()) return;
    if (new_ref != old_ref && live.count(new_ref) != 0) return;
    const O o = it->second;
    live.erase(it);
    if (n > 0) live[new_ref] = {o.l, o.s, p, n, ++clock};
  }
  void quote(TapeRecord& r, Locate l) const {
    std::optional<Price> bp, ap;
    for (const auto& [ref, o] : live) {
      if (o.l != l) continue;
      if (o.s == Side::buy && (!bp || o.p > *bp)) bp = o.p;
      if (o.s == Side::sell && (!ap || o.p < *ap)) ap = o.p;
    }
    std::uint64_t bs = 0, as = 0;
    for (const auto& [ref, o] : live) {
      if (o.l != l) continue;
      if (bp && o.s == Side::buy && o.p == *bp) bs += o.n;
      if (ap && o.s == Side::sell && o.p == *ap) as += o.n;
    }
    constexpr std::uint64_t kMax = std::numeric_limits<std::uint32_t>::max();
    if (bp) {
      r.bid_px = *bp;
      r.bid_sz = static_cast<std::uint32_t>(std::min(bs, kMax));
      if (bs > kMax) r.flags |= kBidSzSaturated;
    }
    if (ap) {
      r.ask_px = *ap;
      r.ask_sz = static_cast<std::uint32_t>(std::min(as, kMax));
      if (as > kMax) r.flags |= kAskSzSaturated;
    }
  }
  bool head(OrderRef r) const {
    const O& o = live.at(r);
    for (const auto& [ref, x] : live) {
      if (x.l == o.l && x.s == o.s && x.p == o.p && x.arrival < o.arrival) return false;
    }
    return true;
  }
  TapeRecord rec(std::uint32_t seq, std::uint64_t ts, char type, std::uint64_t match, Locate l) const {
    TapeRecord r{};
    r.seq = seq;
    r.ts = ts;
    r.type = type;
    r.match = match;
    r.printable = '-';
    r.side = '-';
    r.locate = l;
    quote(r, l);
    r.tstate = tstate[l];
    r.phase = phase;
    if (halted[l]) r.flags |= kOperationalHalt;
    return r;
  }
};

bool same_bytes(const TapeRecord& a, const TapeRecord& b) { return std::memcmp(&a, &b, sizeof a) == 0; }

RC_GTEST_PROP(TapeProperties, FastBookTapeMatchesBruteForceModelAndReferenceBook, ()) {
  Fx<FastBook> fast;
  Fx<ReferenceBook> ref;
  Model m;
  std::vector<TapeRecord> want;
  std::vector<OrderRef> issued;
  OrderRef next = 1;
  std::uint64_t ts = 100;
  const int n = *rc::gen::inRange(1, 400);
  for (int i = 0; i < n; ++i) {
    const auto seq = static_cast<std::uint32_t>(i);
    ts += *rc::gen::inRange<std::uint64_t>(0, 3);  // many ties
    const auto l = *rc::gen::inRange<Locate>(1, 4);
    const auto price = *rc::gen::inRange<Price>(100, 106);
    const auto shares = *rc::gen::inRange<Shares>(1, 200);
    const auto match = *rc::gen::inRange<std::uint64_t>(0, 1'000'000);
    const int kind = *rc::gen::inRange(0, 100);
    const OrderRef target = !issued.empty() && *rc::gen::inRange(0, 20) != 0
                                ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())]
                                : next + 999;
    auto both = [&](auto&& f) {
      f(fast.tb);
      f(ref.tb);
    };
    if (kind < 35) {
      const char sd = *rc::gen::element('B', 'S');
      both([&](auto& tb) { tb.add(l, next, sd, shares, price); });
      m.add(l, next, sd, shares, price);
      issued.push_back(next++);
    } else if (kind < 60) {
      const bool is_c = *rc::gen::inRange(0, 3) == 0;
      const char pr = *rc::gen::element('Y', 'N');
      const auto xp = *rc::gen::inRange<Price>(90, 120);
      const char ty = is_c ? 'C' : 'E';
      const auto it = m.live.find(target);
      TapeRecord r = m.rec(seq, ts, ty, match, it != m.live.end() ? it->second.l : l);
      r.ref = target;
      r.shares = shares;
      r.printable = is_c ? pr : '-';
      if (it != m.live.end()) {
        r.side = it->second.s == Side::buy ? 'B' : 'S';
        r.sign = it->second.s == Side::buy ? std::int8_t{-1} : std::int8_t{1};
        r.at_head = m.head(target) ? 1 : 0;
        r.price = is_c ? xp : it->second.p;
        if (shares > it->second.n) r.flags |= kOverExecute;
      } else {
        r.side = '?';
        r.price = is_c ? xp : 0;
        r.flags |= kOrphan;
      }
      want.push_back(r);
      both([&](auto& tb) { tb.execute(seq, ts, l, target, shares, match, ty, pr, xp); });
      m.reduce(target, shares);
    } else if (kind < 68) {
      both([&](auto& tb) { tb.cancel(l, target, shares); });
      m.reduce(target, shares);
    } else if (kind < 76) {
      both([&](auto& tb) { tb.remove(l, target); });
      m.remove(target);
    } else if (kind < 84) {
      const auto ns = *rc::gen::inRange<Shares>(0, 200);
      const OrderRef nr = *rc::gen::inRange(0, 10) == 0 && !issued.empty()
                              ? issued[*rc::gen::inRange<std::size_t>(0, issued.size())]
                              : next++;
      both([&](auto& tb) { tb.replace(l, target, nr, ns, price); });
      m.replace(target, nr, ns, price);
      issued.push_back(nr);
    } else if (kind < 89) {
      TapeRecord r = m.rec(seq, ts, 'P', match, l);
      r.side = 'B';
      r.shares = shares;
      r.price = price;
      want.push_back(r);
      both([&](auto& tb) { tb.hidden(seq, ts, l, 0, 'B', shares, price, match); });
    } else if (kind < 92) {
      const std::uint64_t qn = *rc::gen::element(std::uint64_t{0}, std::uint64_t{12345}, std::uint64_t{5'000'000'000});
      const char ct = *rc::gen::element('O', 'C', 'H');
      TapeRecord r = m.rec(seq, ts, 'Q', match, l);
      r.side = ct;
      if (qn > std::numeric_limits<std::uint32_t>::max()) {
        r.shares = std::numeric_limits<std::uint32_t>::max();
        r.flags |= kSharesSaturated;
      } else {
        r.shares = static_cast<std::uint32_t>(qn);
      }
      r.price = price;
      want.push_back(r);
      both([&](auto& tb) { tb.cross(seq, ts, l, qn, price, match, ct); });
    } else if (kind < 94) {
      want.push_back(m.rec(seq, ts, 'B', match, l));
      both([&](auto& tb) { tb.broken(seq, ts, l, match); });
    } else if (kind < 96) {
      const char code = *rc::gen::element('O', 'Q', 'M');
      both([&](auto& tb) { tb.system_event(code); });
      if (code == 'Q') m.phase = kRegular;
      if (code == 'M') m.phase = kPost;
    } else if (kind < 98) {
      const char st = *rc::gen::element('T', 'H', 'P', 'Q');
      both([&](auto& tb) { tb.trading_state(l, st); });
      m.tstate[l] = st;
    } else {
      const char mk = *rc::gen::element('Q', 'B', 'X');
      const char act = *rc::gen::element('H', 'T');
      both([&](auto& tb) { tb.operational_halt(l, mk, act); });
      if (mk == 'Q') m.halted[l] = act == 'H';
    }
  }
  RC_ASSERT(fast.recs.size() == want.size());
  RC_ASSERT(ref.recs.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) {
    RC_ASSERT(same_bytes(fast.recs[i], want[i]));
    RC_ASSERT(same_bytes(ref.recs[i], fast.recs[i]));
  }
}

}  // namespace
}  // namespace itchbook::tape
