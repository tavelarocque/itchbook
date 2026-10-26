// exec_probe: read-only measurement of how executions appear in an ITCH day.
// Not part of the repo. Build: c++ -std=c++20 -O2 -I include exec_probe.cpp src/source.cpp src/fast_book.cpp src/reference_book.cpp -lz
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"

using namespace itchbook;
using ull = unsigned long long;

struct Exec {
  std::uint64_t match;
  std::uint64_t ts;
  std::uint32_t shares;
  std::uint16_t locate;
  char type;       // E C P
  char printable;  // C: Y/N, E: '-', P: '-'
  char side;       // E/C: resting side B/S, P: reported side
  std::uint8_t phase;
  std::uint8_t trading;
};

struct Cross {
  std::uint64_t match, ts, shares;
  std::uint16_t locate;
  char type;
  std::uint8_t phase;
};

struct H : itch::NullHandler {
  using NullHandler::on;
  book::FastBook bk;
  std::vector<Exec> ex;
  std::vector<Cross> qs;
  std::uint8_t phase = 0;
  std::array<char, 65536> tstate{};
  std::uint64_t last_ts = 0, ties = 0, regress = 0, msgs = 0;
  std::array<ull, 3> broken{};
  std::vector<std::uint64_t> broken_match;
  // locked/crossed in regular phase while symbol state T
  std::array<std::uint8_t, 65536> cs{};
  std::array<std::uint64_t, 65536> cs_since{};
  ull locked_ep = 0, crossed_ep = 0, one_sided_updates = 0, book_updates_rt = 0;
  double locked_ns = 0, crossed_ns = 0;
  ull locked_lt_1ms = 0, locked_ep_closed = 0;
  // same-timestamp execution runs, regular phase, state T, displayed (E, C printable)
  std::array<std::uint64_t, 65536> run_ts{};
  std::array<char, 65536> run_side{};
  std::array<std::uint32_t, 65536> run_len{};
  std::map<std::uint32_t, ull> run_hist;
  ull mixed_runs = 0, p_same_ts_as_run = 0, p_regular_t = 0;

  H() : bk(std::size_t{1} << 21U) {
    tstate.fill('?');
    ex.reserve(12'000'000);
  }
  template <class M>
  void clock(const M& m) {
    ++msgs;
    const auto ts = m.timestamp();
    if (ts < last_ts) ++regress;
    if (ts == last_ts) ++ties;
    last_ts = ts;
  }
  bool rt(std::uint16_t l) const { return phase == 1 && tstate[l] == 'T'; }
  void book_event(std::uint16_t l, std::uint64_t ts) {
    if (!rt(l)) {
      if (cs[l] != 0) {  // leaving regular/T closes the episode without counting duration
        cs[l] = 0;
      }
      return;
    }
    ++book_updates_rt;
    const auto b = bk.best(l, book::Side::buy);
    const auto a = bk.best(l, book::Side::sell);
    std::uint8_t s = 0;
    if (!b || !a) {
      ++one_sided_updates;
    } else if (b->price > a->price) {
      s = 2;
    } else if (b->price == a->price) {
      s = 1;
    }
    if (s != cs[l]) {
      if (cs[l] == 1) {
        const double d = static_cast<double>(ts - cs_since[l]);
        locked_ns += d;
        ++locked_ep_closed;
        if (d < 1e6) ++locked_lt_1ms;
      } else if (cs[l] == 2) {
        crossed_ns += static_cast<double>(ts - cs_since[l]);
      }
      if (s == 1) ++locked_ep;
      if (s == 2) ++crossed_ep;
      cs[l] = s;
      cs_since[l] = ts;
    }
  }
  void close_run(std::uint16_t l) {
    if (run_len[l] > 0) {
      ++run_hist[std::min<std::uint32_t>(run_len[l], 50)];
      run_len[l] = 0;
    }
  }
  void displayed_exec(std::uint16_t l, std::uint64_t ts, char side) {
    if (!rt(l)) return;
    if (run_len[l] > 0 && run_ts[l] == ts) {
      ++run_len[l];
      if (run_side[l] != side) ++mixed_runs;
    } else {
      close_run(l);
      run_ts[l] = ts;
      run_side[l] = side;
      run_len[l] = 1;
    }
  }

  void on(const itch::SystemEvent& m) {
    clock(m);
    if (m.event_code() == 'Q') { phase = 1; std::printf("system_event\tQ\t%llu\n", (ull)m.timestamp()); }
    if (m.event_code() == 'M') { phase = 2; std::printf("system_event\tM\t%llu\n", (ull)m.timestamp()); }
  }
  void on(const itch::StockTradingAction& m) {
    clock(m);
    tstate[m.stock_locate()] = m.trading_state();
  }
  void on(const itch::AddOrder& m) {
    clock(m);
    bk.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    book_event(m.stock_locate(), m.timestamp());
  }
  void on(const itch::AddOrderMpid& m) {
    clock(m);
    bk.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    book_event(m.stock_locate(), m.timestamp());
  }
  void exec(std::uint64_t ts, std::uint16_t l, std::uint64_t ref, std::uint32_t n, std::uint64_t match, char type,
            char printable) {
    const auto r = bk.resting(ref);
    const char side = r ? (r->side == book::Side::buy ? 'B' : 'S') : '?';
    ex.push_back({match, ts, n, l, type, printable, side, phase, static_cast<std::uint8_t>(tstate[l] == 'T')});
    bk.execute(l, ref, n);
    if (type == 'E' || printable == 'Y') displayed_exec(l, ts, side);
    book_event(l, ts);
  }
  void on(const itch::OrderExecuted& m) {
    clock(m);
    exec(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares(), m.match_number(), 'E', '-');
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    clock(m);
    exec(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares(), m.match_number(), 'C', m.printable());
  }
  void on(const itch::OrderCancel& m) {
    clock(m);
    bk.cancel(m.stock_locate(), m.order_ref(), m.cancelled_shares());
    book_event(m.stock_locate(), m.timestamp());
  }
  void on(const itch::OrderDelete& m) {
    clock(m);
    bk.remove(m.stock_locate(), m.order_ref());
    book_event(m.stock_locate(), m.timestamp());
  }
  void on(const itch::OrderReplace& m) {
    clock(m);
    bk.replace(m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
    book_event(m.stock_locate(), m.timestamp());
  }
  void on(const itch::Trade& m) {
    clock(m);
    const auto l = m.stock_locate();
    ex.push_back({m.match_number(), m.timestamp(), m.shares(), l, 'P', '-', m.side(), phase,
                  static_cast<std::uint8_t>(tstate[l] == 'T')});
    if (rt(l)) {
      ++p_regular_t;
      if (run_len[l] > 0 && run_ts[l] == m.timestamp()) ++p_same_ts_as_run;
    }
  }
  void on(const itch::CrossTrade& m) {
    clock(m);
    qs.push_back({m.match_number(), m.timestamp(), m.shares(), m.stock_locate(), m.cross_type(), phase});
  }
  void on(const itch::BrokenTrade& m) {
    clock(m);
    ++broken[phase];
    broken_match.push_back(m.match_number());
  }
  template <class M>
  void on(const M& m) {
    clock(m);
  }
};

int main(int argc, char** argv) {
  if (argc < 2) return 1;
  auto src = io::open_source(argv[1]);
  auto h = std::make_unique<H>();
  const auto st = itch::parse_stream(*src, *h);
  for (std::size_t l = 0; l < 65536; ++l) h->close_run(static_cast<std::uint16_t>(l));
  const char* ph[3] = {"pre", "regular", "post"};
  std::printf("messages\t%llu\tclean\t%d\n", (ull)st.messages, st.clean() ? 1 : 0);
  std::printf("ts_regressions\t%llu\tts_ties_with_previous_message\t%llu\n", (ull)h->regress, (ull)h->ties);

  std::unordered_set<std::uint64_t> qm;
  for (const auto& q : h->qs) qm.insert(q.match);
  std::unordered_set<std::uint64_t> exm;
  for (const auto& e : h->ex) exm.insert(e.match);

  // executions by type / printable / phase / trading / in-cross / side
  std::map<std::string, std::pair<ull, ull>> tab;  // key -> (count, shares)
  for (const auto& e : h->ex) {
    const bool inq = qm.count(e.match) != 0;
    std::string k = std::string(1, e.type) + "\tprintable=" + e.printable + "\t" + ph[e.phase] + "\tstate_T=" +
                    std::to_string(e.trading) + "\tmatch_in_Q=" + (inq ? "1" : "0") + "\tside=" + e.side;
    auto& v = tab[k];
    ++v.first;
    v.second += e.shares;
  }
  std::printf("\n# executions: type, printable, phase, symbol in state T, match number shared with a Q, side (E/C: resting; P: reported)\n");
  for (const auto& [k, v] : tab) std::printf("exec\t%s\t%llu\t%llu\n", k.c_str(), v.first, v.second);

  std::map<std::string, std::pair<ull, ull>> qtab;
  ull q_with_exec = 0;
  for (const auto& q : h->qs) {
    auto& v = qtab[std::string("type=") + q.type + "\t" + ph[q.phase]];
    ++v.first;
    v.second += q.shares;
    if (exm.count(q.match)) ++q_with_exec;
  }
  std::printf("\n# cross trades Q: cross type, phase, count, shares\n");
  for (const auto& [k, v] : qtab) std::printf("Q\t%s\t%llu\t%llu\n", k.c_str(), v.first, v.second);
  std::printf("Q_sharing_match_with_an_execution\t%llu\tof\t%zu\n", q_with_exec, h->qs.size());

  // match-number uniqueness among executions
  std::unordered_set<std::uint64_t> seen;
  ull dup = 0;
  for (const auto& e : h->ex) {
    if (!seen.insert(e.match).second) ++dup;
  }
  std::printf("\nexecutions\t%zu\tduplicate_match_numbers_among_E_C_P\t%llu\n", h->ex.size(), dup);

  ull bk_e = 0, bk_p = 0, bk_q = 0, bk_none = 0;
  std::unordered_set<std::uint64_t> pm;
  for (const auto& e : h->ex) if (e.type == 'P') pm.insert(e.match);
  for (const auto m : h->broken_match) {
    if (qm.count(m)) ++bk_q;
    else if (pm.count(m)) ++bk_p;
    else if (exm.count(m)) ++bk_e;
    else ++bk_none;
  }
  std::printf("broken_trades_B\tpre\t%llu\tregular\t%llu\tpost\t%llu\trefers_to_E_or_C\t%llu\tP\t%llu\tQ\t%llu\tunmatched\t%llu\n",
              h->broken[0], h->broken[1], h->broken[2], bk_e, bk_p, bk_q, bk_none);

  std::printf("\n# displayed BBO, regular phase, symbol in state T, evaluated after every book event\n");
  std::printf("book_updates\t%llu\tone_sided_updates\t%llu\n", h->book_updates_rt, h->one_sided_updates);
  std::printf("locked_episodes\t%llu\tclosed\t%llu\tunder_1ms\t%llu\tlocked_seconds_total\t%.3f\n", h->locked_ep,
              h->locked_ep_closed, h->locked_lt_1ms, h->locked_ns / 1e9);
  std::printf("crossed_episodes\t%llu\tcrossed_seconds_total\t%.6f\n", h->crossed_ep, h->crossed_ns / 1e9);

  std::printf("\n# runs of displayed executions (E, C printable) on one symbol with the same timestamp, regular, state T\n");
  ull runs = 0, msgs_in = 0;
  for (const auto& [n, c] : h->run_hist) {
    runs += c;
    msgs_in += static_cast<ull>(n) * c;
  }
  std::printf("runs\t%llu\texecutions\t%llu\tmixed_resting_side_runs\t%llu\n", runs, msgs_in, h->mixed_runs);
  for (const auto& [n, c] : h->run_hist) {
    if (n <= 10 || n == 50) std::printf("run_len\t%u%s\t%llu\n", n, n == 50 ? "+" : "", c);
  }
  std::printf("P_regular_T\t%llu\tP_with_same_ts_as_a_displayed_run\t%llu\n", h->p_regular_t, h->p_same_ts_as_run);

  // ---- cross attribution: executions at the same (symbol, timestamp) as a Q print, or within 1 ms before it
  {
    std::unordered_map<std::uint64_t, char> qat;  // key locate<<48 | ts -> cross type
    std::vector<std::vector<std::pair<std::uint64_t, char>>> qby(65536);
    for (const auto& q : h->qs) {
      qat[(std::uint64_t{q.locate} << 48U) | q.ts] = q.type;
      qby[q.locate].push_back({q.ts, q.type});
    }
    std::map<std::string, std::pair<ull, ull>> at, near;
    std::unordered_map<std::uint64_t, ull> mcount;
    for (const auto& e : h->ex) ++mcount[e.match];
    std::map<std::string, ull> dupk;
    for (const auto& e : h->ex) {
      const auto it = qat.find((std::uint64_t{e.locate} << 48U) | e.ts);
      std::string base = std::string(1, e.type) + "\tprintable=" + e.printable;
      char nq = 0;
      for (const auto& [qt, ty] : qby[e.locate]) {
        if (e.ts <= qt && qt - e.ts <= 1'000'000ULL) nq = ty;
      }
      if (it != qat.end()) {
        auto& v = at[base + "\tcross=" + it->second];
        ++v.first; v.second += e.shares;
      } else if (nq != 0) {
        auto& v = near[base + "\tcross=" + nq];
        ++v.first; v.second += e.shares;
      }
      if (mcount[e.match] > 1) {
        ++dupk[base + "\tat_cross_ts=" + (it != qat.end() ? "1" : "0") + "\tphase=" + ph[e.phase]];
      }
    }
    std::printf("\n# executions with the same symbol and timestamp as a Q print (type, printable, cross type O/C/H/I)\n");
    for (const auto& [k, v] : at) std::printf("at_Q_ts\t%s\t%llu\t%llu\n", k.c_str(), v.first, v.second);
    std::printf("# other executions within 1 ms before a Q print on the same symbol\n");
    for (const auto& [k, v] : near) std::printf("within_1ms_before_Q\t%s\t%llu\t%llu\n", k.c_str(), v.first, v.second);
    std::printf("# executions whose match number is shared with another execution\n");
    for (const auto& [k, v] : dupk) std::printf("shared_match\t%s\t%llu\n", k.c_str(), v);
    // Q timing relative to system events
    std::uint64_t qo_min = UINT64_MAX, qo_max = 0, qc_min = UINT64_MAX, qc_max = 0;
    for (const auto& q : h->qs) {
      if (q.type == 'O') { qo_min = std::min(qo_min, q.ts); qo_max = std::max(qo_max, q.ts); }
      if (q.type == 'C') { qc_min = std::min(qc_min, q.ts); qc_max = std::max(qc_max, q.ts); }
    }
    std::printf("Q_open_ts_range_ns\t%llu\t%llu\tQ_close_ts_range_ns\t%llu\t%llu\n", (ull)qo_min, (ull)qo_max, (ull)qc_min, (ull)qc_max);
  }
  return 0;
}
