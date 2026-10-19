#pragma once
// BookBuilder: the ITCH handler that drives a book from a message stream and
// keeps the context needed to judge it: symbol directory, session phase,
// per-symbol trading state, crossed/locked-book monitoring and scheduled
// snapshots.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/reference_book.hpp"
#include "itchbook/itch/messages.hpp"

namespace itchbook::book {

enum class Phase : std::uint8_t { pre_open = 0, regular = 1, post_close = 2 };
inline constexpr std::array<const char*, 3> kPhaseName = {"pre_open", "regular", "post_close"};

struct CrossStats {
  // Index [phase][trading]: trading = 1 if the symbol's last H state was 'T'.
  std::array<std::array<std::uint64_t, 2>, 3> crossed_episodes{};
  std::array<std::array<std::uint64_t, 2>, 3> locked_episodes{};
  std::array<std::array<std::uint64_t, 2>, 3> crossed_updates{};  // book updates that left the book crossed
  struct Example {
    std::uint64_t ts;
    Locate locate;
    Price bid;
    Price ask;
    Phase phase;
    char state;
  };
  std::vector<Example> first_crossed;  // first few crossed episodes, for inspection
  // Crossed updates while the symbol is in state T: how long after the most
  // recent resumption (transition into T) they happened, and on how many symbols.
  std::uint64_t max_crossed_after_resume_ns = 0;
  std::uint64_t symbols_crossed_in_trading = 0;
  std::uint64_t crossed_in_trading_never_halted = 0;  // crossed in T on a symbol with no resumption: real anomaly
};

struct ExecStats {
  // Executions of orders that were not first in their queue, by message type
  // (0 = E, 1 = C) and phase; and all executions for the denominators.
  std::array<std::array<std::uint64_t, 3>, 2> behind_head{};
  std::array<std::array<std::uint64_t, 3>, 2> total{};
};

struct Snapshot {
  std::uint64_t at_ns;  // state after every message with timestamp <= at_ns
  std::uint64_t checksum;
  std::uint64_t live_orders;
  std::uint64_t levels;
  struct Top {
    Locate locate;
    Side side;
    std::vector<LevelInfo> levels;
  };
  std::vector<Top> tops;
};

struct BuilderOptions {
  std::vector<std::uint64_t> snapshot_times_ns;  // ascending
  std::vector<std::string> top_symbols;          // symbols whose top levels each snapshot records
  std::size_t top_depth = 5;
};

template <class BookT = ReferenceBook>
class BookBuilder {
 public:
  using Options = BuilderOptions;

  explicit BookBuilder(Options opts) : opts_(std::move(opts)) { trading_.fill('?'); }

  template <class M>
  void on(const M& m) {
    advance_clock(m.timestamp());
    handle(m);
  }

  // Call after the stream ends: emits any snapshots not yet reached.
  void finish() {
    while (next_snap_ < opts_.snapshot_times_ns.size()) {
      take_snapshot(opts_.snapshot_times_ns[next_snap_++]);
    }
  }

  [[nodiscard]] const BookT& book() const noexcept { return book_; }
  [[nodiscard]] const CrossStats& cross_stats() const noexcept { return cross_; }
  [[nodiscard]] const ExecStats& exec_stats() const noexcept { return exec_; }
  [[nodiscard]] const std::vector<Snapshot>& snapshots() const noexcept { return snaps_; }
  [[nodiscard]] std::string_view symbol(Locate l) const noexcept { return symbols_[l]; }
  [[nodiscard]] std::size_t symbol_count() const noexcept { return n_symbols_; }

 private:
  void advance_clock(std::uint64_t ts) {
    while (next_snap_ < opts_.snapshot_times_ns.size() && ts > opts_.snapshot_times_ns[next_snap_]) {
      take_snapshot(opts_.snapshot_times_ns[next_snap_++]);
    }
  }

  void take_snapshot(std::uint64_t at) {
    Snapshot s{at, book_.checksum(), book_.live_orders(), book_.level_count(), {}};
    for (const Locate l : top_locates_) {
      for (const Side side : {Side::buy, Side::sell}) {
        s.tops.push_back({l, side, book_.levels(l, side, opts_.top_depth)});
      }
    }
    snaps_.push_back(std::move(s));
  }

  void check_cross(Locate l, std::uint64_t ts) {
    const auto bid = book_.best(l, Side::buy);
    const auto ask = book_.best(l, Side::sell);
    const int state = (bid && ask) ? (bid->price > ask->price ? 2 : (bid->price == ask->price ? 1 : 0)) : 0;
    const auto p = static_cast<std::size_t>(phase_);
    const std::size_t t = trading_[l] == 'T' ? 1 : 0;
    if (state == 2) {
      ++cross_.crossed_updates[p][t];
      if (t == 1) {
        if (!crossed_in_trading_[l]) {
          crossed_in_trading_[l] = true;
          ++cross_.symbols_crossed_in_trading;
        }
        if (resumed_at_[l] == kNever) {
          ++cross_.crossed_in_trading_never_halted;
        } else {
          cross_.max_crossed_after_resume_ns = std::max(cross_.max_crossed_after_resume_ns, ts - resumed_at_[l]);
        }
      }
    }
    if (state != cross_state_[l]) {
      if (state == 2) {
        ++cross_.crossed_episodes[p][t];
        if (cross_.first_crossed.size() < 25 && bid && ask) {
          cross_.first_crossed.push_back({ts, l, bid->price, ask->price, phase_, trading_[l]});
        }
      } else if (state == 1) {
        ++cross_.locked_episodes[p][t];
      }
      cross_state_[l] = static_cast<std::uint8_t>(state);
    }
  }

  // ---- message handlers -----------------------------------------------------
  void handle(const itch::SystemEvent& m) {
    if (m.event_code() == 'Q') {
      phase_ = Phase::regular;
    } else if (m.event_code() == 'M') {
      phase_ = Phase::post_close;
    }
  }
  void handle(const itch::StockDirectory& m) {
    std::string_view s = m.stock();
    while (!s.empty() && s.back() == ' ') {
      s.remove_suffix(1);
    }
    if (symbols_[m.stock_locate()].empty()) {
      ++n_symbols_;
    }
    symbols_[m.stock_locate()] = std::string(s);
    for (const auto& want : opts_.top_symbols) {
      if (want == s) {
        top_locates_.push_back(m.stock_locate());
      }
    }
  }
  void handle(const itch::StockTradingAction& m) {
    const Locate l = m.stock_locate();
    // A resumption is a move into T from a state seen earlier (not the
    // start-of-day spin, where the previous state is unknown).
    if (m.trading_state() == 'T' && trading_[l] != 'T' && trading_[l] != '?') {
      resumed_at_[l] = m.timestamp();
    }
    trading_[l] = m.trading_state();
  }
  void count_exec(std::size_t type, OrderRef ref) {
    const auto p = static_cast<std::size_t>(phase_);
    ++exec_.total[type][p];
    if (!book_.at_head(ref)) {
      ++exec_.behind_head[type][p];
    }
  }
  void handle(const itch::AddOrder& m) {
    book_.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::AddOrderMpid& m) {
    book_.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::OrderExecuted& m) {
    count_exec(0, m.order_ref());
    book_.execute(m.stock_locate(), m.order_ref(), m.executed_shares());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::OrderExecutedWithPrice& m) {
    count_exec(1, m.order_ref());
    book_.execute(m.stock_locate(), m.order_ref(), m.executed_shares());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::OrderCancel& m) {
    book_.cancel(m.stock_locate(), m.order_ref(), m.cancelled_shares());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::OrderDelete& m) {
    book_.remove(m.stock_locate(), m.order_ref());
    check_cross(m.stock_locate(), m.timestamp());
  }
  void handle(const itch::OrderReplace& m) {
    book_.replace(m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
    check_cross(m.stock_locate(), m.timestamp());
  }
  template <class M>
  void handle(const M& /*m*/) {}  // P, Q, B, I, L, ...: no effect on the displayed book

  Options opts_;
  BookT book_;
  CrossStats cross_;
  ExecStats exec_;
  static constexpr std::uint64_t kNever = UINT64_MAX;
  std::vector<std::uint64_t> resumed_at_ = std::vector<std::uint64_t>(65536, kNever);
  std::vector<bool> crossed_in_trading_ = std::vector<bool>(65536, false);
  std::vector<Snapshot> snaps_;
  std::size_t next_snap_ = 0;
  Phase phase_ = Phase::pre_open;
  std::array<char, 65536> trading_{};
  std::array<std::uint8_t, 65536> cross_state_{};  // 0 normal, 1 locked, 2 crossed
  std::vector<std::string> symbols_ = std::vector<std::string>(65536);
  std::size_t n_symbols_ = 0;
  std::vector<Locate> top_locates_;
};

}  // namespace itchbook::book
