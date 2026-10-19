// libFuzzer target: arbitrary bytes through the frame splitter and every
// field accessor. The first input byte picks the source chunk size so the
// fuzzer also explores frames split across reads. Invariants are checked with
// __builtin_trap so a violation is reported as a crash.

#include <cstddef>
#include <cstdint>
#include <span>

#include "itchbook/io/source.hpp"
#include "itchbook/itch/fields.hpp"
#include "itchbook/itch/stream.hpp"

namespace {

struct TouchAll : itchbook::itch::NullHandler {
  using NullHandler::on;
  std::uint64_t acc = 0;
  std::uint64_t n = 0;
  template <class M>
  void on(const M& m) {
    ++n;
    itchbook::itch::visit_fields(m, [&](std::string_view, std::uint64_t v) { acc += v; });
  }
};

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
  if (size == 0) {
    return 0;
  }
  const std::size_t chunk = std::size_t{data[0]} + 1;
  const std::span<const std::uint8_t> body(data + 1, size - 1);
  itchbook::io::MemorySource src(body, chunk);
  TouchAll h;
  const auto st = itchbook::itch::parse_stream(src, h, {}, 64);
  if (st.bytes + st.truncated_tail_bytes != body.size()) {
    __builtin_trap();
  }
  if (st.messages != h.n || st.frames != st.messages + st.unknown_type + st.bad_length + st.end_of_session_frames) {
    __builtin_trap();
  }
  return 0;
}
