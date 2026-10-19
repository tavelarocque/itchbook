// Field coverage: for every message type, setting any single byte (other
// than the type byte) must change exactly one visited field. This proves the
// accessors tile each message with no gaps and no overlaps, which is what
// makes the per-field fingerprint a complete check of the decoder.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include "itchbook/itch/fields.hpp"
#include "itchbook/itch/messages.hpp"

namespace itchbook::itch {
namespace {

template <class M>
void check_tiling() {
  for (std::size_t pos = 1; pos < M::kSize; ++pos) {
    std::vector<std::uint8_t> b(M::kSize, 0);
    b[0] = static_cast<std::uint8_t>(M::kType);
    b[pos] = 0x01;
    int nonzero = 0;
    std::string which;
    visit_fields(M{b.data()}, [&](std::string_view name, std::uint64_t v) {
      if (v != 0) {
        ++nonzero;
        which += std::string(name) + " ";
      }
    });
    EXPECT_EQ(nonzero, 1) << "type " << M::kType << " byte " << pos << " hits: " << which;
  }
}

template <class... Ms>
void check_all(TypeList<Ms...> /*unused*/) {
  (check_tiling<Ms>(), ...);
}

TEST(Fields, EveryByteBelongsToExactlyOneField) { check_all(AllMessages{}); }

TEST(Fields, AlphaEncodesBigEndian) {
  EXPECT_EQ(alpha_as_u64("AAPL    "), 0x4141504C20202020ULL);
  EXPECT_EQ(alpha_as_u64("N"), 0x4EULL);
}

}  // namespace
}  // namespace itchbook::itch
