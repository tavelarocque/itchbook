// itch-book: rebuild the full L3 book for a trading day with ReferenceBook and
// report its health.
//
// stdout is canonical, sorted text that tools/itchref/bookref.py reproduces
// independently (counters, a whole-market checksum at each snapshot time, the
// top levels of selected symbols); `diff` compares them. stderr carries the
// diagnostics that have no Python counterpart: crossed/locked-book episodes by
// session phase, end-of-day residue and timing.
//
//   itch-book [--engine fast|reference] [--every MINUTES] [--symbols AAPL,MSFT] [--depth N] [--max-messages N] FILE|-

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "itchbook/book/builder.hpp"
#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)
using ull = unsigned long long;

std::string hms(std::uint64_t ns) {
  const std::uint64_t s = ns / 1'000'000'000ULL;
  char out[16];
  std::snprintf(out, sizeof out, "%02llu:%02llu:%02llu", static_cast<ull>(s / 3600), static_cast<ull>((s / 60) % 60),
                static_cast<ull>(s % 60));
  return out;
}

std::vector<std::string> split(std::string_view s, char sep) {
  std::vector<std::string> out;
  while (!s.empty()) {
    const auto k = s.find(sep);
    out.emplace_back(s.substr(0, k));
    if (k == std::string_view::npos) {
      break;
    }
    s.remove_prefix(k + 1);
  }
  return out;
}

void print_stats(const book::BookStats& s) {
  const std::pair<const char*, std::uint64_t> rows[] = {
      {"adds", s.adds},
      {"executes", s.executes},
      {"cancels", s.cancels},
      {"deletes", s.deletes},
      {"replaces", s.replaces},
      {"orphan_execute", s.orphan_execute},
      {"orphan_cancel", s.orphan_cancel},
      {"orphan_delete", s.orphan_delete},
      {"orphan_replace", s.orphan_replace},
      {"duplicate_add", s.duplicate_add},
      {"duplicate_replace", s.duplicate_replace},
      {"bad_side", s.bad_side},
      {"over_execute", s.over_execute},
      {"over_cancel", s.over_cancel},
      {"locate_mismatch", s.locate_mismatch},
      {"replace_to_zero", s.replace_to_zero},
      {"executed_at_head", s.executed_at_head},
      {"executed_behind_head", s.executed_behind_head},
      {"cancel_to_zero", s.cancel_to_zero},
      {"max_live_orders", s.max_live_orders},
  };
  for (const auto& [k, v] : rows) {
    std::printf("stat\t%s\t%llu\n", k, static_cast<ull>(v));
  }
}

template <class BookT>
int run(const std::string& path, book::BuilderOptions opts, const itch::StreamLimits& limits) {
  auto src = io::open_source(path);
  book::BookBuilder<BookT> b(std::move(opts));
  const auto t0 = std::chrono::steady_clock::now();
  const auto st = itch::parse_stream(*src, b, limits);
  b.finish();
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const auto& bk = b.book();

  std::printf("messages\t%llu\n", static_cast<ull>(st.messages));
  print_stats(bk.stats());
  std::printf("end_live_orders\t%llu\n", static_cast<ull>(bk.live_orders()));
  std::printf("end_levels\t%llu\n", static_cast<ull>(bk.level_count()));
  for (const auto& s : b.snapshots()) {
    std::printf("snap\t%s\t%llu\t%llu\t%llu\n", hms(s.at_ns).c_str(), static_cast<ull>(s.checksum),
                static_cast<ull>(s.live_orders), static_cast<ull>(s.levels));
    auto tops = s.tops;
    std::sort(tops.begin(), tops.end(), [&](const auto& x, const auto& y) {
      return std::pair(b.symbol(x.locate), x.side) < std::pair(b.symbol(y.locate), y.side);
    });
    for (const auto& t : tops) {
      for (std::size_t r = 0; r < t.levels.size(); ++r) {
        const auto& l = t.levels[r];
        std::printf("top\t%s\t%.*s\t%c\t%zu\t%u\t%llu\t%u\n", hms(s.at_ns).c_str(),
                    static_cast<int>(b.symbol(t.locate).size()), b.symbol(t.locate).data(), book::side_char(t.side), r,
                    l.price, static_cast<ull>(l.shares), l.orders);
      }
    }
  }

  // ---- diagnostics (stderr) --------------------------------------------------
  const auto& s = bk.stats();
  std::fprintf(stderr, "itch-book (%s): %s\n  %.1f s, %.0f ns/msg (includes decompression), %zu symbols\n",
               std::is_same_v<BookT, book::FastBook> ? "fast" : "reference", path.c_str(), secs,
               secs * 1e9 / static_cast<double>(st.messages != 0 ? st.messages : 1), b.symbol_count());
  std::fprintf(stderr, "  orphans %llu, anomalies %llu, max live orders %llu, end live %zu, end levels %zu\n",
               static_cast<ull>(s.orphans()), static_cast<ull>(s.anomalies()), static_cast<ull>(s.max_live_orders),
               bk.live_orders(), bk.level_count());
  const auto ex = static_cast<double>(s.executed_at_head + s.executed_behind_head);
  std::fprintf(stderr, "  executions at head of queue: %.4f%% of %.0f\n",
               ex > 0 ? 100.0 * static_cast<double>(s.executed_at_head) / ex : 0.0, ex);
  const auto& c = b.cross_stats();
  std::fprintf(stderr, "  crossed/locked episodes by phase (symbol state T | other):\n");
  for (std::size_t p = 0; p < 3; ++p) {
    std::fprintf(stderr, "    %-10s crossed %llu | %llu   locked %llu | %llu   crossed updates %llu | %llu\n",
                 book::kPhaseName[p], static_cast<ull>(c.crossed_episodes[p][1]),
                 static_cast<ull>(c.crossed_episodes[p][0]), static_cast<ull>(c.locked_episodes[p][1]),
                 static_cast<ull>(c.locked_episodes[p][0]), static_cast<ull>(c.crossed_updates[p][1]),
                 static_cast<ull>(c.crossed_updates[p][0]));
  }
  std::fprintf(stderr,
               "    crossed while in state T: %llu symbols, max %.3f s after a resumption, %llu updates on symbols "
               "never resumed\n",
               static_cast<ull>(c.symbols_crossed_in_trading), static_cast<double>(c.max_crossed_after_resume_ns) / 1e9,
               static_cast<ull>(c.crossed_in_trading_never_halted));
  const auto& x = b.exec_stats();
  for (std::size_t ty = 0; ty < 2; ++ty) {
    std::fprintf(stderr, "  %s executions behind head of queue:", ty == 0 ? "E" : "C");
    for (std::size_t p = 0; p < 3; ++p) {
      std::fprintf(stderr, "  %s %llu/%llu", book::kPhaseName[p], static_cast<ull>(x.behind_head[ty][p]),
                   static_cast<ull>(x.total[ty][p]));
    }
    std::fprintf(stderr, "\n");
  }
  for (const auto& e : c.first_crossed) {
    const auto sym = b.symbol(e.locate);
    std::fprintf(stderr, "    crossed: %s %.*s bid %u ask %u phase %s state %c\n", hms(e.ts).c_str(),
                 static_cast<int>(sym.size()), sym.data(), e.bid, e.ask,
                 book::kPhaseName[static_cast<std::size_t>(e.phase)], e.state);
  }
  if (!st.clean() && !st.stopped_at_limit) {
    std::fprintf(stderr, "  WARNING: stream not clean\n");
    return 2;
  }
  return 0;
}

void usage() {
  std::fprintf(stderr,
               "usage: itch-book [--engine fast|reference] [--every MINUTES] [--symbols A,B] [--depth N] "
               "[--max-messages N] FILE|-\n");
}

template <class T>
bool parse_num(std::string_view v, T& out) {
  return std::from_chars(v.data(), v.data() + v.size(), out).ec == std::errc{};
}

}  // namespace

int main(int argc, char** argv) {
  book::BuilderOptions opts;
  bool fast = true;
  opts.top_symbols = {"AAPL", "MSFT", "SPY", "QQQ", "TSLA", "AMD", "ZVZZT"};
  itch::StreamLimits limits;
  std::uint64_t every_min = 30;
  std::string path;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    const bool has_value = i + 1 < argc;
    if (a == "--every" && has_value) {
      if (!parse_num(argv[++i], every_min) || every_min == 0) {
        usage();
        return 1;
      }
    } else if (a == "--engine" && has_value) {
      const std::string_view e = argv[++i];
      if (e != "fast" && e != "reference") {
        usage();
        return 1;
      }
      fast = e == "fast";
    } else if (a == "--symbols" && has_value) {
      opts.top_symbols = split(argv[++i], ',');
    } else if (a == "--depth" && has_value) {
      if (!parse_num(argv[++i], opts.top_depth)) {
        usage();
        return 1;
      }
    } else if (a == "--max-messages" && has_value) {
      if (!parse_num(argv[++i], limits.max_messages)) {
        usage();
        return 1;
      }
    } else if (path.empty() && (a == "-" || !a.starts_with("-"))) {
      path = a;
    } else {
      usage();
      return 1;
    }
  }
  if (path.empty()) {
    usage();
    return 1;
  }
  // Snapshots on a fixed wall-clock grid from 03:30 to 20:30 ET.
  constexpr std::uint64_t kMin = 60'000'000'000ULL;
  for (std::uint64_t t = 210 * kMin; t <= 1230 * kMin; t += every_min * kMin) {
    opts.snapshot_times_ns.push_back(t);
  }
  try {
    return fast ? run<book::FastBook>(path, std::move(opts), limits)
                : run<book::ReferenceBook>(path, std::move(opts), limits);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-book: %s\n", e.what());
    return 1;
  }
}
