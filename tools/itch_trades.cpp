// itch-trades: the per-execution tape of a day.
//
//   itch-trades [--engine fast|reference] [--max-messages N] --out PREFIX FILE|-
//
// Replays the full L3 book and writes, in file order, one 64-byte record per
// E, C, P, Q and B message (layout: include/itchbook/tape/record.hpp), with the
// displayed quote read before each message is applied. Nothing is filtered.
//
// Outputs:
//   PREFIX.tape.gz       gzip: FileHeader + records
//   PREFIX.symbols.tsv   the stock directory (locates are per day)
//   PREFIX.meta.json     counts, stream health, engine
#include <bit>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <zlib.h>

#include "itchbook/book/fast_book.hpp"
#include "itchbook/book/reference_book.hpp"
#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"
#include "itchbook/tape/builder.hpp"
#include "itchbook/tape/record.hpp"

using namespace itchbook;
using ull = unsigned long long;

static_assert(std::endian::native == std::endian::little, "tape records are written in host order");

namespace {

class GzWriter {
 public:
  explicit GzWriter(const std::string& path) : f_(gzopen(path.c_str(), "wb1")) {
    if (f_ == nullptr) throw std::runtime_error("cannot open " + path);
    buf_.reserve(kBuf);
    tape::FileHeader h{};
    for (int i = 0; i < 8; ++i) h.magic[i] = tape::kMagic[i];
    h.version = tape::kVersion;
    h.record_size = sizeof(tape::TapeRecord);
    h.endian = 0x01020304U;
    write(&h, sizeof h);
  }
  GzWriter(const GzWriter&) = delete;
  GzWriter& operator=(const GzWriter&) = delete;
  ~GzWriter() {
    if (f_ != nullptr) gzclose(f_);
  }
  void push(const tape::TapeRecord& r) {
    buf_.push_back(r);
    if (buf_.size() == kBuf) flush();
  }
  void close() {
    flush();
    if (gzclose(f_) != Z_OK) throw std::runtime_error("gzclose failed");
    f_ = nullptr;
  }

 private:
  static constexpr std::size_t kBuf = 1U << 16U;
  void flush() {
    if (!buf_.empty()) write(buf_.data(), buf_.size() * sizeof(tape::TapeRecord));
    buf_.clear();
  }
  void write(const void* p, std::size_t n) {
    if (gzwrite(f_, p, static_cast<unsigned>(n)) != static_cast<int>(n)) throw std::runtime_error("gzwrite failed");
  }
  gzFile f_;
  std::vector<tape::TapeRecord> buf_;
};

struct Counts {
  ull records = 0, E = 0, C = 0, P = 0, Q = 0, B = 0;
  ull orphan = 0, over_execute = 0, saturated = 0;
};

struct Sink {
  GzWriter* w;
  Counts* c;
  void operator()(const tape::TapeRecord& r) const {
    w->push(r);
    ++c->records;
    switch (r.type) {
      case 'E': ++c->E; break;
      case 'C': ++c->C; break;
      case 'P': ++c->P; break;
      case 'Q': ++c->Q; break;
      case 'B': ++c->B; break;
      default: break;
    }
    if ((r.flags & tape::kOrphan) != 0) ++c->orphan;
    if ((r.flags & tape::kOverExecute) != 0) ++c->over_execute;
    if ((r.flags & (tape::kSharesSaturated | tape::kBidSzSaturated | tape::kAskSzSaturated)) != 0) ++c->saturated;
  }
};

struct DirEntry {
  std::string symbol;
  char market_category, issue_classification, etp_flag;
};

template <class Book>
struct Handler : itch::NullHandler {
  using NullHandler::on;
  tape::TapeBuilder<Book, Sink> tb;
  std::vector<std::pair<std::uint16_t, DirEntry>> dir;
  std::uint64_t seq = 0;
  bool seq_overflow = false;

  Handler(Book& b, Sink s) : tb(b, s) {}

  std::uint32_t next() {
    if (seq > UINT32_MAX) seq_overflow = true;
    return static_cast<std::uint32_t>(seq++);
  }
  void on(const itch::SystemEvent& m) {
    next();
    tb.system_event(m.event_code());
  }
  void on(const itch::StockDirectory& m) {
    next();
    auto s = m.stock();
    while (!s.empty() && s.back() == ' ') s.remove_suffix(1);
    dir.push_back({m.stock_locate(), {std::string(s), m.market_category(), m.issue_classification(), m.etp_flag()}});
  }
  void on(const itch::StockTradingAction& m) {
    next();
    tb.trading_state(m.stock_locate(), m.trading_state());
  }
  void on(const itch::OperationalHalt& m) {
    next();
    tb.operational_halt(m.stock_locate(), m.market_code(), m.halt_action());
  }
  void on(const itch::AddOrder& m) {
    next();
    tb.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::AddOrderMpid& m) {
    next();
    tb.add(m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::OrderExecuted& m) {
    const auto s = next();
    tb.execute(s, m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares(), m.match_number(), 'E', '-', 0);
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    const auto s = next();
    tb.execute(s, m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares(), m.match_number(), 'C',
               m.printable(), m.execution_price());
  }
  void on(const itch::OrderCancel& m) {
    next();
    tb.cancel(m.stock_locate(), m.order_ref(), m.cancelled_shares());
  }
  void on(const itch::OrderDelete& m) {
    next();
    tb.remove(m.stock_locate(), m.order_ref());
  }
  void on(const itch::OrderReplace& m) {
    next();
    tb.replace(m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
  }
  void on(const itch::Trade& m) {
    const auto s = next();
    tb.hidden(s, m.timestamp(), m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price(), m.match_number());
  }
  void on(const itch::CrossTrade& m) {
    const auto s = next();
    tb.cross(s, m.timestamp(), m.stock_locate(), m.shares(), m.cross_price(), m.match_number(), m.cross_type());
  }
  void on(const itch::BrokenTrade& m) {
    const auto s = next();
    tb.broken(s, m.timestamp(), m.stock_locate(), m.match_number());
  }
  template <class M>
  void on(const M&) {
    next();
  }
};

template <class Book>
int run(const std::string& engine, const std::string& in, const std::string& prefix, const itch::StreamLimits& lim) {
  auto book = std::make_unique<Book>();
  Counts c;
  GzWriter w(prefix + ".tape.gz");
  auto h = std::make_unique<Handler<Book>>(*book, Sink{&w, &c});
  auto src = io::open_source(in);
  const auto st = itch::parse_stream(*src, *h, lim);
  w.close();

  if (std::FILE* f = std::fopen((prefix + ".symbols.tsv").c_str(), "w")) {
    std::fprintf(f, "locate\tsymbol\tmarket_category\tissue_classification\tetp_flag\n");
    for (const auto& [l, d] : h->dir) {
      std::fprintf(f, "%u\t%s\t%c\t%c\t%c\n", static_cast<unsigned>(l), d.symbol.c_str(), d.market_category,
                   d.issue_classification, d.etp_flag);
    }
    std::fclose(f);
  } else {
    throw std::runtime_error("cannot write symbols.tsv");
  }
  const auto& bs = book->stats();
  if (std::FILE* f = std::fopen((prefix + ".meta.json").c_str(), "w")) {
    std::fprintf(f,
                 "{\n  \"source\": \"%s\",\n  \"engine\": \"%s\",\n  \"tape_version\": %u,\n  \"record_size\": %zu,\n"
                 "  \"messages\": %llu,\n  \"clean\": %s,\n  \"stopped_at_limit\": %s,\n  \"seq_overflow\": %s,\n"
                 "  \"records\": %llu,\n  \"E\": %llu,\n  \"C\": %llu,\n  \"P\": %llu,\n  \"Q\": %llu,\n  \"B\": %llu,\n"
                 "  \"orphan_executions\": %llu,\n  \"over_executions\": %llu,\n  \"saturated_records\": %llu,\n"
                 "  \"book_anomalies\": %llu,\n  \"symbols\": %zu\n}\n",
                 in.c_str(), engine.c_str(), tape::kVersion, sizeof(tape::TapeRecord), static_cast<ull>(st.messages),
                 st.clean() ? "true" : "false", st.stopped_at_limit ? "true" : "false",
                 h->seq_overflow ? "true" : "false", c.records, c.E, c.C, c.P, c.Q, c.B, c.orphan, c.over_execute,
                 c.saturated, static_cast<ull>(bs.anomalies()), h->dir.size());
    std::fclose(f);
  } else {
    throw std::runtime_error("cannot write meta.json");
  }
  std::fprintf(stderr, "itch-trades (%s): %llu messages, %llu records (E %llu, C %llu, P %llu, Q %llu, B %llu)%s\n",
               engine.c_str(), static_cast<ull>(st.messages), c.records, c.E, c.C, c.P, c.Q, c.B,
               st.clean() || st.stopped_at_limit ? "" : "  WARNING: stream not clean");
  return st.clean() || st.stopped_at_limit ? 0 : 1;
}

bool parse_num(std::string_view s, std::uint64_t& out) {
  const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
  return r.ec == std::errc{} && r.ptr == s.data() + s.size();
}

int usage() {
  std::fprintf(stderr, "usage: itch-trades [--engine fast|reference] [--max-messages N] --out PREFIX FILE|-\n");
  return 2;
}

}  // namespace

int main(int argc, char** argv) {
  std::string engine = "fast", out, in;
  itch::StreamLimits lim;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    const bool has_value = i + 1 < argc;
    if (a == "--engine" && has_value) {
      engine = argv[++i];
      if (engine != "fast" && engine != "reference") return usage();
    } else if (a == "--max-messages" && has_value) {
      if (!parse_num(argv[++i], lim.max_messages)) return usage();
    } else if (a == "--out" && has_value) {
      out = argv[++i];
    } else if (in.empty() && (a == "-" || !a.starts_with("--"))) {
      in = a;
    } else {
      return usage();
    }
  }
  if (in.empty() || out.empty()) return usage();
  try {
    return engine == "fast" ? run<book::FastBook>(engine, in, out, lim)
                            : run<book::ReferenceBook>(engine, in, out, lim);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-trades: %s\n", e.what());
    return 1;
  }
}
