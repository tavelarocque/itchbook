// itch-lob: fixed-grid order-book features for the M5 forecasting study.
//
//   itch-lob --rank [--top N] FILE
//       Activity per symbol during regular hours (trade prints, shares, book
//       events, VWAP), most active first. Used once, on the first day, to fix
//       the study's universe.
//
//   itch-lob --symbols A,B,C | --universe FILE  [--grid-ms 250] --out PREFIX FILE
//       Writes PREFIX.SYMBOL.npy (int32, rows x 56, layout in
//       include/itchbook/lob/recorder.hpp) for each symbol, and PREFIX.meta.json.

#include <algorithm>
#include <bit>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"
#include "itchbook/lob/recorder.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)
using ull = unsigned long long;
using book::Locate;
using book::OrderRef;
using book::Price;
using book::Shares;

std::string trim(std::string_view s) {
  while (!s.empty() && s.back() == ' ') {
    s.remove_suffix(1);
  }
  return std::string(s);
}

std::vector<std::string> split(std::string_view s, char sep) {
  std::vector<std::string> out;
  while (!s.empty()) {
    const auto k = s.find(sep);
    if (k != 0) {
      out.emplace_back(s.substr(0, k));
    }
    if (k == std::string_view::npos) {
      break;
    }
    s.remove_prefix(k + 1);
  }
  return out;
}

// ---- recording --------------------------------------------------------------------

struct Recorder : itch::NullHandler {
  using NullHandler::on;
  lob::LobRecorder rec;
  explicit Recorder(lob::RecorderOptions o) : rec(std::move(o)) {}

  template <class M>
  void tick(const M& m) {
    rec.clock(m.timestamp());
  }
  void on(const itch::StockDirectory& m) {
    tick(m);
    rec.directory(m.stock_locate(), trim(m.stock()));
  }
  void on(const itch::StockTradingAction& m) {
    tick(m);
    rec.trading_state(m.stock_locate(), m.trading_state());
  }
  void on(const itch::AddOrder& m) {
    tick(m);
    rec.add(m.timestamp(), m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::AddOrderMpid& m) {
    tick(m);
    rec.add(m.timestamp(), m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::OrderExecuted& m) {
    tick(m);
    rec.execute(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares());
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    tick(m);
    rec.execute(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares(), m.printable() == 'Y');
  }
  void on(const itch::OrderCancel& m) {
    tick(m);
    rec.cancel(m.timestamp(), m.stock_locate(), m.order_ref(), m.cancelled_shares());
  }
  void on(const itch::OrderDelete& m) {
    tick(m);
    rec.remove(m.timestamp(), m.stock_locate(), m.order_ref());
  }
  void on(const itch::OrderReplace& m) {
    tick(m);
    rec.replace(m.timestamp(), m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
  }
  void on(const itch::Trade& m) {
    tick(m);
    rec.hidden_trade(m.timestamp(), m.stock_locate(), m.side(), m.shares());
  }
  template <class M>
  void on(const M& m) {
    tick(m);
  }
};

void write_npy(const std::string& path, const std::vector<lob::Row>& rows) {
  std::string header = "{'descr': '<i4', 'fortran_order': False, 'shape': (" + std::to_string(rows.size()) + ", " +
                       std::to_string(lob::kCols) + "), }";
  // Magic (6) + version (2) + header length (2) + header, padded with spaces to a multiple of 64, ending in \n.
  const std::size_t total = ((10 + header.size() + 1 + 63) / 64) * 64;
  header.append(total - 10 - header.size() - 1, ' ');
  header += '\n';
  std::ofstream f(path, std::ios::binary);
  if (!f) {
    throw std::runtime_error("cannot write " + path);
  }
  const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0};
  f.write(magic, sizeof magic);
  const auto hl = static_cast<std::uint16_t>(header.size());
  const char len[2] = {static_cast<char>(hl & 0xFFU), static_cast<char>(hl >> 8U)};
  f.write(len, 2);
  f.write(header.data(), static_cast<std::streamsize>(header.size()));
  static_assert(std::endian::native == std::endian::little, "npy writer assumes a little-endian host");
  f.write(reinterpret_cast<const char*>(rows.data()),  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
          static_cast<std::streamsize>(rows.size() * sizeof(lob::Row)));
  if (!f) {
    throw std::runtime_error("short write to " + path);
  }
}

// ---- ranking ----------------------------------------------------------------------

struct Ranker : itch::NullHandler {
  using NullHandler::on;
  book::FastBook book;
  struct Stat {
    std::uint64_t prints = 0;
    std::uint64_t shares = 0;
    double notional = 0;
    std::uint64_t events = 0;
  };
  std::vector<Stat> st = std::vector<Stat>(65536);
  std::vector<std::string> sym = std::vector<std::string>(65536);
  std::vector<char> etp = std::vector<char>(65536, ' ');

  static bool regular(std::uint64_t ts) { return ts > lob::kOpenNs && ts <= lob::kCloseNs; }
  void on(const itch::StockDirectory& m) {
    sym[m.stock_locate()] = trim(m.stock());
    etp[m.stock_locate()] = m.etp_flag();
  }
  void count(std::uint64_t ts, Locate l) {
    if (regular(ts)) {
      ++st[l].events;
    }
  }
  void on(const itch::AddOrder& m) {
    book.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    count(m.timestamp(), m.stock_locate());
  }
  void on(const itch::AddOrderMpid& m) {
    book.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
    count(m.timestamp(), m.stock_locate());
  }
  void exec(std::uint64_t ts, OrderRef ref, Shares n, std::optional<Price> px) {
    const auto r = book.resting(ref);
    if (r && regular(ts)) {
      Stat& s = st[r->locate];
      const auto v = std::min(n, r->shares);
      ++s.prints;
      s.shares += v;
      s.notional += static_cast<double>(v) * static_cast<double>(px.value_or(r->price)) / 1e4;
      ++s.events;
    }
  }
  void on(const itch::OrderExecuted& m) {
    exec(m.timestamp(), m.order_ref(), m.executed_shares(), std::nullopt);
    book.execute(m.stock_locate(), m.order_ref(), m.executed_shares());
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    if (m.printable() == 'Y') {
      exec(m.timestamp(), m.order_ref(), m.executed_shares(), m.execution_price());
    }
    book.execute(m.stock_locate(), m.order_ref(), m.executed_shares());
  }
  void on(const itch::OrderCancel& m) {
    book.cancel(m.stock_locate(), m.order_ref(), m.cancelled_shares());
    count(m.timestamp(), m.stock_locate());
  }
  void on(const itch::OrderDelete& m) {
    book.remove(m.stock_locate(), m.order_ref());
    count(m.timestamp(), m.stock_locate());
  }
  void on(const itch::OrderReplace& m) {
    book.replace(m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
    count(m.timestamp(), m.stock_locate());
  }
};

void usage() {
  std::fprintf(stderr,
               "usage: itch-lob --rank [--top N] FILE\n"
               "       itch-lob (--symbols A,B | --universe FILE) [--grid-ms MS] --out PREFIX FILE\n");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    bool rank = false;
    std::size_t top = 100;
    std::vector<std::string> symbols;
    std::uint64_t grid_ms = 250;
    std::string out;
    std::string path;
    for (int i = 1; i < argc; ++i) {
      const std::string_view a = argv[i];
      const bool v = i + 1 < argc;
      if (a == "--rank") {
        rank = true;
      } else if (a == "--top" && v) {
        top = std::stoul(argv[++i]);
      } else if (a == "--symbols" && v) {
        symbols = split(argv[++i], ',');
      } else if (a == "--universe" && v) {
        std::ifstream f(argv[++i]);
        std::string line;
        while (std::getline(f, line)) {
          if (!line.empty() && line[0] != '#') {
            symbols.push_back(split(line, '\t').at(0));
          }
        }
      } else if (a == "--grid-ms" && v) {
        grid_ms = std::stoull(argv[++i]);
      } else if (a == "--out" && v) {
        out = argv[++i];
      } else if (path.empty() && (a == "-" || !a.starts_with("-"))) {
        path = a;
      } else {
        usage();
        return 1;
      }
    }
    if (path.empty() || (!rank && (symbols.empty() || out.empty() || grid_ms == 0))) {
      usage();
      return 1;
    }
    const auto t0 = std::chrono::steady_clock::now();
    auto src = io::open_source(path);
    if (rank) {
      auto h = std::make_unique<Ranker>();
      itch::parse_stream(*src, *h);
      std::vector<std::size_t> order;
      for (std::size_t l = 0; l < h->st.size(); ++l) {
        if (h->st[l].prints > 0) {
          order.push_back(l);
        }
      }
      std::sort(order.begin(), order.end(), [&](auto a, auto b) { return h->st[a].prints > h->st[b].prints; });
      std::printf("symbol\tetp\tprints\tshares\tnotional_usd\tvwap\tbook_events\n");
      for (std::size_t i = 0; i < std::min(top, order.size()); ++i) {
        const auto& s = h->st[order[i]];
        std::printf("%s\t%c\t%llu\t%llu\t%.0f\t%.4f\t%llu\n", h->sym[order[i]].c_str(), h->etp[order[i]],
                    static_cast<ull>(s.prints), static_cast<ull>(s.shares), s.notional,
                    s.shares > 0 ? s.notional / static_cast<double>(s.shares) : 0.0, static_cast<ull>(s.events));
      }
      return 0;
    }
    lob::RecorderOptions o;
    o.symbols = symbols;
    o.grid_ns = grid_ms * 1'000'000ULL;
    auto h = std::make_unique<Recorder>(o);
    const itch::StreamStats st = itch::parse_stream(*src, *h);
    h->rec.finish();
    std::size_t missing = 0;
    std::string meta = "{\n  \"source\": \"" + path + "\",\n  \"grid_ns\": " + std::to_string(o.grid_ns) +
                       ",\n  \"start_ns\": " + std::to_string(o.start_ns) +
                       ",\n  \"rows\": " + std::to_string(h->rec.rows()) +
                       ",\n  \"depth\": " + std::to_string(lob::kDepth) +
                       ",\n  \"messages\": " + std::to_string(st.messages) + ",\n  \"symbols\": [";
    for (std::size_t i = 0; i < symbols.size(); ++i) {
      if (!h->rec.seen(i)) {
        ++missing;
        std::fprintf(stderr, "itch-lob: %s not in this day's directory\n", symbols[i].c_str());
      } else {
        write_npy(out + "." + symbols[i] + ".npy", h->rec.rows_of(i));
      }
      meta += (i > 0 ? ", \"" : "\"") + symbols[i] + "\"";
    }
    meta += "],\n  \"missing\": " + std::to_string(missing) + "\n}\n";
    std::ofstream(out + ".meta.json") << meta;
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("messages\t%llu\nsymbols\t%zu\nmissing\t%zu\nrows\t%llu\nseconds\t%.1f\n",
                static_cast<ull>(st.messages), symbols.size(), missing, static_cast<ull>(h->rec.rows()), secs);
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-lob: %s\n", e.what());
    return 1;
  }
}
