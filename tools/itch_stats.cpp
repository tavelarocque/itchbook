// itch-stats: stream an ITCH 5.0 file and print per-type message counts and,
// optionally, a per-field fingerprint (the sum mod 2^64 of every field of
// every message, grouped by type).
//
// The output on stdout is a canonical, sorted, tab-separated text that the
// independent Python reference (tools/itchref) prints in the same format, so
// the two can be compared with `diff`. Timing goes to stderr.
//
//   itch-stats [--fingerprint] [--max-messages N] FILE|-

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/fields.hpp"
#include "itchbook/itch/stream.hpp"

namespace {

using namespace itchbook;  // NOLINT(google-build-using-namespace)

struct TypeAcc {
  std::uint64_t count = 0;
  std::vector<std::pair<std::string_view, std::uint64_t>> sums;
};

template <bool kFingerprint>
struct StatsHandler : itch::NullHandler {
  std::array<TypeAcc, 256> acc{};
  std::uint64_t first_ts = UINT64_MAX;
  std::uint64_t last_ts = 0;
  std::uint64_t ts_regressions = 0;  // timestamp lower than the previous message's
  std::vector<std::pair<char, std::uint64_t>> system_events;
  std::uint64_t max_order_ref = 0;  // largest reference introduced by A, F or U
  std::uint16_t max_locate = 0;

  template <class M>
  void on(const M& m) {
    TypeAcc& a = acc[static_cast<std::uint8_t>(M::kType)];
    ++a.count;
    const std::uint64_t ts = m.timestamp();
    if (first_ts == UINT64_MAX) {
      first_ts = ts;
    }
    if (ts < last_ts) {
      ++ts_regressions;
    }
    last_ts = ts;
    if constexpr (std::is_same_v<M, itch::SystemEvent>) {
      system_events.emplace_back(m.event_code(), ts);
    }
    if constexpr (std::is_same_v<M, itch::AddOrder> || std::is_same_v<M, itch::AddOrderMpid>) {
      max_order_ref = std::max(max_order_ref, m.order_ref());
    }
    if constexpr (std::is_same_v<M, itch::OrderReplace>) {
      max_order_ref = std::max(max_order_ref, m.new_order_ref());
    }
    max_locate = std::max(max_locate, m.stock_locate());
    if constexpr (kFingerprint) {
      std::size_t i = 0;
      itch::visit_fields(m, [&](std::string_view name, std::uint64_t v) {
        if (i == a.sums.size()) {
          a.sums.emplace_back(name, 0);
        }
        a.sums[i++].second += v;  // wraps mod 2^64 by definition
      });
    }
  }
};

std::string hms(std::uint64_t ns) {
  const std::uint64_t s = ns / 1'000'000'000ULL;
  char out[32];
  std::snprintf(out, sizeof out, "%02llu:%02llu:%02llu.%09llu", static_cast<unsigned long long>(s / 3600),
                static_cast<unsigned long long>((s / 60) % 60), static_cast<unsigned long long>(s % 60),
                static_cast<unsigned long long>(ns % 1'000'000'000ULL));
  return out;
}

template <bool kFingerprint>
int run(const std::string& path, const itch::StreamLimits& limits) {
  auto src = io::open_source(path);
  StatsHandler<kFingerprint> h;
  const auto t0 = std::chrono::steady_clock::now();
  const itch::StreamStats st = itch::parse_stream(*src, h, limits);
  const auto t1 = std::chrono::steady_clock::now();
  const double secs = std::chrono::duration<double>(t1 - t0).count();

  std::printf("frames\t%llu\n", static_cast<unsigned long long>(st.frames));
  std::printf("messages\t%llu\n", static_cast<unsigned long long>(st.messages));
  std::printf("end_of_session_frames\t%llu\n", static_cast<unsigned long long>(st.end_of_session_frames));
  std::printf("unknown_type\t%llu\n", static_cast<unsigned long long>(st.unknown_type));
  std::printf("bad_length\t%llu\n", static_cast<unsigned long long>(st.bad_length));
  std::printf("truncated_tail_bytes\t%llu\n", static_cast<unsigned long long>(st.truncated_tail_bytes));
  for (std::size_t t = 0; t < h.acc.size(); ++t) {
    if (h.acc[t].count != 0) {
      std::printf("count\t%c\t%llu\n", static_cast<char>(t), static_cast<unsigned long long>(h.acc[t].count));
    }
  }
  if constexpr (kFingerprint) {
    for (std::size_t t = 0; t < h.acc.size(); ++t) {
      for (const auto& [name, sum] : h.acc[t].sums) {
        std::printf("fp\t%c\t%.*s\t%llu\n", static_cast<char>(t), static_cast<int>(name.size()), name.data(),
                    static_cast<unsigned long long>(sum));
      }
    }
  }

  std::fprintf(stderr, "itch-stats: %s\n", path.c_str());
  std::fprintf(stderr, "  %.3f s, %.1f M msg/s, %.0f MB/s raw, %.1f ns/msg (includes decompression)\n", secs,
               static_cast<double>(st.messages) / secs / 1e6, static_cast<double>(st.bytes) / secs / 1e6,
               secs * 1e9 / static_cast<double>(st.messages == 0 ? 1 : st.messages));
  if (h.first_ts != UINT64_MAX) {
    std::fprintf(stderr, "  timestamps %s .. %s, %llu regressions\n", hms(h.first_ts).c_str(), hms(h.last_ts).c_str(),
                 static_cast<unsigned long long>(h.ts_regressions));
  }
  std::fprintf(stderr, "  max order reference %llu, max stock locate %u\n",
               static_cast<unsigned long long>(h.max_order_ref), static_cast<unsigned>(h.max_locate));
  std::fprintf(stderr, "  system events:");
  for (std::size_t i = 0; i < h.system_events.size() && i < 16; ++i) {
    const auto& [code, ts] = h.system_events[i];
    const auto u = static_cast<unsigned char>(code);
    if (u >= 0x20 && u < 0x7f) {
      std::fprintf(stderr, " %c@%s", code, hms(ts).c_str());
    } else {
      std::fprintf(stderr, " \\x%02x@%s", u, hms(ts).c_str());
    }
  }
  std::fprintf(stderr, h.system_events.size() > 16 ? " ...\n" : "\n");
  if (st.stopped_at_limit) {
    std::fprintf(stderr, "  stopped at --max-messages\n");
  }
  if (st.source_truncated) {
    std::fprintf(stderr, "  WARNING: compressed input is truncated (partial file)\n");
  }
  if (!st.clean()) {
    std::fprintf(stderr, "  WARNING: stream not clean (see counters above)\n");
    return st.stopped_at_limit && st.unknown_type == 0 && st.bad_length == 0 ? 0 : 2;
  }
  return 0;
}

void usage() { std::fprintf(stderr, "usage: itch-stats [--fingerprint] [--max-messages N] FILE|-\n"); }

}  // namespace

int main(int argc, char** argv) {
  bool fingerprint = false;
  itch::StreamLimits limits;
  std::string path;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    if (a == "--fingerprint") {
      fingerprint = true;
    } else if (a == "--max-messages" && i + 1 < argc) {
      const std::string_view v = argv[++i];
      if (std::from_chars(v.data(), v.data() + v.size(), limits.max_messages).ec != std::errc{}) {
        usage();
        return 1;
      }
    } else if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else if (path.empty()) {
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
  try {
    return fingerprint ? run<true>(path, limits) : run<false>(path, limits);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "itch-stats: %s\n", e.what());
    return 1;
  }
}
