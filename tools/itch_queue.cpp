// itch-queue: the M4 queue-position study on a trading day.
//
// Follows a sample of orders that join the best bid or ask behind other
// orders during regular hours and, at every trade print at their price,
// compares each L2 model's estimate of the shares ahead with the truth from
// the L3 book. Writes, under the --out prefix:
//
//   PREFIX.orders.tsv.gz   one row per tracked order: entry state, real
//                          outcome, and per model the estimation error and
//                          the fill it would have predicted
//   PREFIX.points.tsv.gz   for a sub-sample of orders, every evaluation point
//                          (truth and each model's estimate before a print)
//   PREFIX.cancels.tsv     where cancellations come from, by the tracked
//                          order's true relative position; also the
//                          calibration input for --calibration
//
//   itch-queue [--sample N] [--events N] [--verify N] [--calibration FILE]
//              [--max-messages N] --out PREFIX FILE|-

#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <zlib.h>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/stream.hpp"
#include "itchbook/queue/study.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)
using ull = unsigned long long;

class GzOut {
 public:
  explicit GzOut(const std::string& path) : f_(gzopen(path.c_str(), "wb1")) {
    if (f_ == nullptr) {
      throw std::runtime_error("cannot write " + path);
    }
    gzbuffer(f_, 1U << 20U);
  }
  GzOut(const GzOut&) = delete;
  GzOut& operator=(const GzOut&) = delete;
  ~GzOut() {
    if (f_ != nullptr) {
      gzclose(f_);
    }
  }
  void write(std::string_view s) {
    if (!s.empty() && gzwrite(f_, s.data(), static_cast<unsigned>(s.size())) == 0) {
      throw std::runtime_error("gzwrite failed");
    }
  }

 private:
  gzFile f_;
};

void num(std::string& out, std::uint64_t v) {
  char b[24];
  const auto r = std::to_chars(b, b + sizeof b, v);
  out.append(b, r.ptr);
}
void num(std::string& out, std::int64_t v) {
  char b[24];
  const auto r = std::to_chars(b, b + sizeof b, v);
  out.append(b, r.ptr);
}
void num(std::string& out, double v) {
  char b[32];
  const int n = std::snprintf(b, sizeof b, "%.6g", v);
  out.append(b, static_cast<std::size_t>(n));
}

struct Handler : itch::NullHandler {
  using NullHandler::on;
  queue::QueueStudy study;
  std::vector<std::string> symbols = std::vector<std::string>(65536);

  explicit Handler(queue::StudyOptions o) : study(std::move(o)) {}

  void on(const itch::SystemEvent& m) {
    if (m.event_code() == 'Q') {
      study.set_active(true, m.timestamp());
    } else if (m.event_code() == 'M') {
      study.set_active(false, m.timestamp());
    }
  }
  void on(const itch::StockDirectory& m) {
    std::string_view s = m.stock();
    while (!s.empty() && s.back() == ' ') {
      s.remove_suffix(1);
    }
    symbols[m.stock_locate()] = std::string(s);
  }
  void on(const itch::StockTradingAction& m) { study.set_trading(m.stock_locate(), m.trading_state() == 'T'); }
  void on(const itch::AddOrder& m) {
    study.add(m.timestamp(), m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::AddOrderMpid& m) {
    study.add(m.timestamp(), m.stock_locate(), m.order_ref(), m.side(), m.shares(), m.price());
  }
  void on(const itch::OrderExecuted& m) {
    study.execute(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares());
  }
  void on(const itch::OrderExecutedWithPrice& m) {
    study.execute(m.timestamp(), m.stock_locate(), m.order_ref(), m.executed_shares());
  }
  void on(const itch::OrderCancel& m) {
    study.cancel(m.timestamp(), m.stock_locate(), m.order_ref(), m.cancelled_shares());
  }
  void on(const itch::OrderDelete& m) { study.remove(m.timestamp(), m.stock_locate(), m.order_ref()); }
  void on(const itch::OrderReplace& m) {
    study.replace(m.timestamp(), m.stock_locate(), m.original_order_ref(), m.new_order_ref(), m.shares(), m.price());
  }
};

queue::Calibration load_calibration(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read " + path);
  }
  std::vector<double> phi;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.starts_with("bin")) {
      continue;
    }
    std::istringstream ss(line);
    std::string bin;
    std::string lo;
    std::string hi;
    std::string events;
    std::string total;
    std::string ahead;
    double p = 0;
    ss >> bin >> lo >> hi >> events >> total >> ahead >> p;
    phi.push_back(p);
  }
  if (phi.size() != queue::Calibration::kBins) {
    throw std::runtime_error(path + ": expected " + std::to_string(queue::Calibration::kBins) + " bins");
  }
  return queue::Calibration(phi);
}

void usage() {
  std::fprintf(stderr,
               "usage: itch-queue [--sample N] [--events N] [--verify N] [--calibration FILE] "
               "[--max-messages N] --out PREFIX FILE|-\n");
}

bool parse_u64(std::string_view v, std::uint64_t& out) {
  return std::from_chars(v.data(), v.data() + v.size(), out).ec == std::errc{};
}

}  // namespace

int main(int argc, char** argv) {
  try {
    queue::StudyOptions o;
    std::string path;
    std::string out;
    std::string cal_path;
    itch::StreamLimits limits;
    for (int i = 1; i < argc; ++i) {
      const std::string_view a = argv[i];
      const bool v = i + 1 < argc;
      if (a == "--sample" && v) {
        if (!parse_u64(argv[++i], o.sample)) {
          usage();
          return 1;
        }
      } else if (a == "--events" && v) {
        if (!parse_u64(argv[++i], o.event_sample)) {
          usage();
          return 1;
        }
      } else if (a == "--verify" && v) {
        if (!parse_u64(argv[++i], o.verify_every)) {
          usage();
          return 1;
        }
      } else if (a == "--max-messages" && v) {
        if (!parse_u64(argv[++i], limits.max_messages)) {
          usage();
          return 1;
        }
      } else if (a == "--calibration" && v) {
        cal_path = argv[++i];
      } else if (a == "--out" && v) {
        out = argv[++i];
      } else if (path.empty() && (a == "-" || !a.starts_with("-"))) {
        path = a;
      } else {
        usage();
        return 1;
      }
    }
    if (path.empty() || out.empty()) {
      usage();
      return 1;
    }
    if (!cal_path.empty()) {
      o.calibration = load_calibration(cal_path);
    }

    GzOut orders(out + ".orders.tsv.gz");
    GzOut points(out + ".points.tsv.gz");
    std::string buf;
    buf.reserve(1U << 16U);
    buf =
        "ref\tlocate\tside\tprice\tadd_ts\tend_ts\tend\town0\tahead0\torders_ahead0\tspread\tfirst_fill_ts"
        "\tfill_trades\tfilled\ttruth_at_fill\tn_eval\tn_trades";
    for (const auto name : queue::kModelName) {
      for (const char* f : {"abs", "err", "err_fill", "pred_ts", "pred_trades", "pred_filled"}) {
        buf += '\t';
        buf += name;
        buf += '_';
        buf += f;
      }
    }
    buf += '\n';
    orders.write(buf);
    buf = "ref\tts\tothers\ttruth";
    for (const auto name : queue::kModelName) {
      buf += '\t';
      buf += name;
    }
    buf += '\n';
    points.write(buf);

    std::uint64_t n_records = 0;
    std::uint64_t n_points = 0;
    std::string row;
    o.on_record = [&](const queue::OrderRecord& r) {
      row.clear();
      num(row, std::uint64_t{r.ref});
      row += '\t';
      num(row, std::uint64_t{r.locate});
      row += r.side == book::Side::buy ? "\tB\t" : "\tS\t";
      num(row, std::uint64_t{r.price});
      row += '\t';
      num(row, r.add_ts);
      row += '\t';
      num(row, r.end_ts);
      row += '\t';
      row += queue::kEndName[static_cast<std::size_t>(r.end)];
      for (const std::uint64_t x : {std::uint64_t{r.own0}, r.ahead0, std::uint64_t{r.orders_ahead0},
                                    std::uint64_t{r.spread}, r.first_fill_ts, r.fill_trades, r.filled}) {
        row += '\t';
        num(row, x);
      }
      row += '\t';
      num(row, r.truth_at_fill);
      row += '\t';
      num(row, r.n_eval);
      row += '\t';
      num(row, r.n_trades);
      for (const queue::ModelResult& m : r.m) {
        row += '\t';
        num(row, m.sum_abs_err);
        row += '\t';
        num(row, m.sum_err);
        row += '\t';
        num(row, m.err_at_fill);
        row += '\t';
        num(row, m.pred_fill_ts);
        row += '\t';
        num(row, m.pred_fill_trades);
        row += '\t';
        num(row, m.pred_filled);
      }
      row += '\n';
      buf += row;
      if (buf.size() > (1U << 16U)) {
        orders.write(buf);
        buf.clear();
      }
      ++n_records;
    };
    buf.clear();
    std::string pbuf;
    o.on_point = [&](const queue::EvalPoint& p) {
      num(pbuf, std::uint64_t{p.ref});
      pbuf += '\t';
      num(pbuf, p.ts);
      pbuf += '\t';
      num(pbuf, p.others);
      pbuf += '\t';
      num(pbuf, p.truth);
      for (const double e : p.est) {
        pbuf += '\t';
        num(pbuf, e);
      }
      pbuf += '\n';
      if (pbuf.size() > (1U << 16U)) {
        points.write(pbuf);
        pbuf.clear();
      }
      ++n_points;
    };
    const bool calibrated = o.calibration.fitted();
    auto h = std::make_unique<Handler>(std::move(o));

    const auto t0 = std::chrono::steady_clock::now();
    auto src = io::open_source(path);
    const itch::StreamStats st = itch::parse_stream(*src, *h, limits);
    h->study.set_active(false, 0);  // flush anything still tracked (truncated input)
    orders.write(buf);
    points.write(pbuf);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const queue::CancelLocation& cl = h->study.cancel_location();
    if (FILE* f = std::fopen((out + ".cancels.tsv").c_str(), "w")) {
      std::fprintf(
          f, "# share of cancelled volume that was ahead of a tracked order, by its true position r = ahead/others\n");
      std::fprintf(f, "bin\tlo\thi\tevents\tvolume\tvolume_ahead\tphi\n");
      for (std::size_t b = 0; b < queue::Calibration::kBins; ++b) {
        const double lo = static_cast<double>(b) / queue::Calibration::kBins;
        const double phi = cl.total[b] > 0 ? cl.ahead[b] / cl.total[b] : queue::Calibration::centre(b);
        std::fprintf(f, "%zu\t%.2f\t%.2f\t%llu\t%.0f\t%.0f\t%.6f\n", b, lo, lo + 1.0 / queue::Calibration::kBins,
                     static_cast<ull>(cl.events[b]), cl.total[b], cl.ahead[b], phi);
      }
      std::fclose(f);
    }
    const queue::QueueFlow& qf = h->study.queue_flow();
    if (FILE* f = std::fopen((out + ".flow.tsv").c_str(), "w")) {
      std::fprintf(f, "# shares that left the queue around tracked orders, regular hours, by cause\n");
      std::fprintf(f, "cause\tahead\tbehind\n");
      for (std::size_t c = 0; c < queue::kCauses; ++c) {
        std::fprintf(f, "%s\t%.0f\t%.0f\n", queue::kCauseName[c], qf.ahead[c], qf.behind[c]);
      }
      std::fprintf(f, "replace_same_price_ahead_net_visible\t%.0f\t\n", qf.replace_same_net_visible);
      std::fclose(f);
    }
    std::printf("messages\t%llu\n", static_cast<ull>(st.messages));
    std::printf("tracked_orders\t%llu\n", static_cast<ull>(n_records));
    std::printf("eval_points_written\t%llu\n", static_cast<ull>(n_points));
    std::printf("calibration\t%s\n", calibrated ? cal_path.c_str() : "none (calibrated model = uniform)");
    std::printf("verified\t%llu\tfailures\t%llu\n", static_cast<ull>(h->study.verified()),
                static_cast<ull>(h->study.verify_failures()));
    const book::BookStats& bs = h->study.book().stats();
    const std::uint64_t orphans = bs.orphan_execute + bs.orphan_cancel + bs.orphan_delete + bs.orphan_replace;
    std::printf("orphans\t%llu\n", static_cast<ull>(orphans));
    std::printf("seconds\t%.1f\n", secs);
    return h->study.verify_failures() == 0 ? 0 : 2;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-queue: %s\n", e.what());
    return 1;
  }
}
