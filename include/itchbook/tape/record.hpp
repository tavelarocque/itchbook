#pragma once
// The per-execution tape: one fixed-size record per execution-related ITCH
// message (E, C, P, Q, B), in file order. Written by itch-trades, read by
// tools/tape/tape.py (whose numpy dtype must match this layout).
//
// File layout: a 24-byte FileHeader, then N TapeRecords, the whole file gzip
// compressed. Integers are little-endian (the writer asserts a little-endian
// host). Nothing is filtered: auction fills, halts, hidden prints and broken
// trades are all recorded, with the fields needed to filter them downstream.

#include <cstddef>
#include <cstdint>

namespace itchbook::tape {

inline constexpr char kMagic[8] = {'I', 'T', 'C', 'H', 'T', 'A', 'P', 'E'};
inline constexpr std::uint32_t kVersion = 1;

struct FileHeader {
  char magic[8];
  std::uint32_t version;
  std::uint32_t record_size;
  std::uint32_t endian;  // 0x01020304 as written by the host
  std::uint32_t reserved;
};
static_assert(sizeof(FileHeader) == 24);

// Session phase from the S system events: before 'Q', between 'Q' and 'M', after 'M'.
inline constexpr std::uint8_t kPre = 0;
inline constexpr std::uint8_t kRegular = 1;
inline constexpr std::uint8_t kPost = 2;

// flags bits
inline constexpr std::uint8_t kSharesSaturated = 1U << 0U;  // shares did not fit in 32 bits (Q only)
inline constexpr std::uint8_t kBidSzSaturated = 1U << 1U;
inline constexpr std::uint8_t kAskSzSaturated = 1U << 2U;
inline constexpr std::uint8_t kOverExecute = 1U << 3U;      // E/C for more shares than were resting
inline constexpr std::uint8_t kOrphan = 1U << 4U;           // E/C against an unknown order reference
inline constexpr std::uint8_t kOperationalHalt = 1U << 5U;  // NASDAQ 'h' halt in effect for the symbol

struct TapeRecord {
  std::uint64_t ts;      // ns since midnight, US Eastern
  std::uint64_t match;   // match number (all types)
  std::uint64_t ref;     // order reference (E, C, P as sent; 0 for Q, B)
  std::uint32_t seq;     // 0-based index of the message among decoded messages: file order
  std::uint32_t shares;  // shares in the message (E/C: as sent, not clamped; Q: saturated)
  std::uint32_t price;   // E: resting order's price; C: execution price; P: message price;
                         // Q: cross price; B: 0. ITCH Price(4) units.
  std::uint32_t bid_px;  // displayed best bid BEFORE the message is applied; 0 = side empty
  std::uint32_t ask_px;  // displayed best ask BEFORE the message; 0 = side empty
  std::uint32_t bid_sz;  // shares at the best bid (saturated at 2^32 - 1)
  std::uint32_t ask_sz;  // shares at the best ask (saturated)
  std::uint16_t locate;  // E/C: the resting order's locate (header locate if orphan); else header
  char type;             // 'E' 'C' 'P' 'Q' 'B'
  char printable;        // C: 'Y'/'N' as sent; E, P, Q, B: '-'
  char side;             // E/C: resting side 'B'/'S' ('?' if orphan); P: side as reported (always
                         // 'B' since 14 Jul 2014); Q: cross type 'O'/'C'/'H'/'I'; B: '-'
  std::int8_t sign;      // E/C: +1 buyer-initiated (resting sell), -1 seller-initiated; else 0
  char tstate;           // symbol's last H trading state at the message ('?' if none yet)
  std::uint8_t phase;    // kPre, kRegular, kPost
  std::uint8_t at_head;  // E/C: 1 if the resting order was first in its queue before the fill
  std::uint8_t flags;    // bit set above
  std::uint8_t pad[2];
};
static_assert(sizeof(TapeRecord) == 64);
static_assert(offsetof(TapeRecord, seq) == 24);
static_assert(offsetof(TapeRecord, locate) == 52);
static_assert(offsetof(TapeRecord, flags) == 61);

}  // namespace itchbook::tape
