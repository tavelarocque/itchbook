#include "itchbook/book/fast_book.hpp"

#include <bit>

#include <limits>
#include "itchbook/itch/endian.hpp"

namespace itchbook::book {

FastBook::FastBook(std::size_t expected_orders)
    : books_(std::size_t{std::numeric_limits<Locate>::max()} + 1), index_(expected_orders) {
  slots_.reserve(expected_orders);
}

// ---- levels -----------------------------------------------------------------

namespace {
// Index of the first level that is not worse than p, scanning from the back
// (best end) first because almost all activity is near the top of the book.
template <class Vec, class Better>
std::size_t lower(const Vec& v, Price p, Better better) {
  std::size_t i = v.size();
  for (int steps = 0; steps < 8 && i > 0; ++steps) {
    if (!better(v[i - 1].price, p)) {
      // v[i-1] is not better than p: answer is i-1 if it equals p's rank, else i.
      return better(p, v[i - 1].price) ? i : i - 1;
    }
    --i;
  }
  // Long book: binary search over [0, i) for the first element not worse than p.
  std::size_t lo = 0;
  std::size_t hi = i;
  while (lo < hi) {
    const std::size_t mid = lo + (hi - lo) / 2;
    if (better(p, v[mid].price)) {
      lo = mid + 1;  // v[mid] is worse than p
    } else {
      hi = mid;
    }
  }
  return lo;
}
}  // namespace

std::size_t FastBook::position(const Side_& v, Side s, Price p) const noexcept {
  return lower(v, p, [s](Price a, Price b) { return better(s, a, b); });
}

std::uint32_t FastBook::level_for_insert(Locate l, Side s, Price p) {
  Side_& v = books_[l].sides[static_cast<int>(s)];
  const std::size_t i = position(v, s, p);
  if (i < v.size() && v[i].price == p) {
    return v[i].level;
  }
  std::uint32_t li = free_level_;
  if (li != kNil) {
    free_level_ = levels_[li].head;
  } else {
    levels_.emplace_back();
    li = static_cast<std::uint32_t>(levels_.size() - 1);
  }
  levels_[li] = Level{0, p, 0, kNil, kNil};
  v.insert(v.begin() + static_cast<std::ptrdiff_t>(i), PriceRef{p, li});
  ++level_count_;
  return li;
}

void FastBook::remove_level(Locate l, Side s, Price p, std::uint32_t level) {
  Side_& v = books_[l].sides[static_cast<int>(s)];
  v.erase(v.begin() + static_cast<std::ptrdiff_t>(position(v, s, p)));
  levels_[level].head = free_level_;
  free_level_ = level;
  --level_count_;
}

// ---- orders -----------------------------------------------------------------

std::uint32_t FastBook::alloc_slot() {
  if (free_ != kNil) {
    const std::uint32_t s = free_;
    free_ = slots_[s].next;
    return s;
  }
  slots_.emplace_back();
  return static_cast<std::uint32_t>(slots_.size() - 1);
}

void FastBook::insert(OrderRef ref, Locate locate, Side side, Shares shares, Price price) {
  const std::uint32_t s = alloc_slot();
  const std::uint32_t li = level_for_insert(locate, side, price);
  Level& lvl = levels_[li];
  slots_[s] = Slot{ref, price, shares, lvl.tail, kNil, li, locate, side};
  if (lvl.tail != kNil) {
    slots_[lvl.tail].next = s;
  } else {
    lvl.head = s;
  }
  lvl.tail = s;
  ++lvl.count;
  lvl.shares += shares;
  index_.insert(ref, s);
  ++live_;
  stats_.max_live_orders = std::max<std::uint64_t>(stats_.max_live_orders, live_);
}

void FastBook::erase(std::uint32_t s) {
  const Slot o = slots_[s];
  Level* lvl = &levels_[o.level];
  if (o.prev != kNil) {
    slots_[o.prev].next = o.next;
  } else {
    lvl->head = o.next;
  }
  if (o.next != kNil) {
    slots_[o.next].prev = o.prev;
  } else {
    lvl->tail = o.prev;
  }
  lvl->shares -= o.shares;
  if (--lvl->count == 0) {
    remove_level(o.locate, o.side, o.price, o.level);
  }
  index_.erase(o.ref);
  slots_[s].next = free_;
  free_ = s;
  --live_;
}

void FastBook::reduce(std::uint32_t s, Shares n) {
  Slot& o = slots_[s];
  if (n >= o.shares) {
    erase(s);
    return;
  }
  levels_[o.level].shares -= n;
  o.shares -= n;
}

void FastBook::add(Locate locate, OrderRef ref, char side, Shares shares, Price price) {
  if (side != 'B' && side != 'S') {
    ++stats_.bad_side;
    return;
  }
  if (find(ref) != kNil) {
    ++stats_.duplicate_add;
    return;
  }
  ++stats_.adds;
  insert(ref, locate, side == 'B' ? Side::buy : Side::sell, shares, price);
}

void FastBook::execute(Locate locate, OrderRef ref, Shares shares) {
  const std::uint32_t s = find(ref);
  if (s == kNil) {
    ++stats_.orphan_execute;
    return;
  }
  ++stats_.executes;
  const Slot& o = slots_[s];
  if (locate != o.locate) {
    ++stats_.locate_mismatch;
  }
  ++(o.prev == kNil ? stats_.executed_at_head : stats_.executed_behind_head);
  if (shares > o.shares) {
    ++stats_.over_execute;
  }
  reduce(s, shares);
}

void FastBook::cancel(Locate locate, OrderRef ref, Shares shares) {
  const std::uint32_t s = find(ref);
  if (s == kNil) {
    ++stats_.orphan_cancel;
    return;
  }
  ++stats_.cancels;
  const Slot& o = slots_[s];
  if (locate != o.locate) {
    ++stats_.locate_mismatch;
  }
  if (shares > o.shares) {
    ++stats_.over_cancel;
  }
  if (shares >= o.shares) {
    ++stats_.cancel_to_zero;
  }
  reduce(s, shares);
}

void FastBook::remove(Locate locate, OrderRef ref) {
  const std::uint32_t s = find(ref);
  if (s == kNil) {
    ++stats_.orphan_delete;
    return;
  }
  ++stats_.deletes;
  if (locate != slots_[s].locate) {
    ++stats_.locate_mismatch;
  }
  erase(s);
}

void FastBook::replace(Locate locate, OrderRef old_ref, OrderRef new_ref, Shares shares, Price price) {
  const std::uint32_t s = find(old_ref);
  if (s == kNil) {
    ++stats_.orphan_replace;
    return;
  }
  if (new_ref != old_ref && find(new_ref) != kNil) {
    ++stats_.duplicate_replace;
    return;
  }
  ++stats_.replaces;
  const Locate l = slots_[s].locate;
  const Side side = slots_[s].side;
  if (locate != l) {
    ++stats_.locate_mismatch;
  }
  erase(s);
  if (shares == 0) {
    ++stats_.replace_to_zero;
    return;
  }
  insert(new_ref, l, side, shares, price);
}

// ---- prefetch ---------------------------------------------------------------

void FastBook::prefetch(const std::uint8_t* p, std::size_t len, int stage) const noexcept {
  const char type = static_cast<char>(p[0]);
  const bool add = type == 'A' || type == 'F';
  const bool by_ref = type == 'E' || type == 'C' || type == 'X' || type == 'D' || type == 'U';
  if ((!add && !by_ref) || len < 20) {
    return;
  }
  const OrderRef ref = itch::load_be64(p + 11);
  if (stage == 0) {
    index_.prefetch(ref);
    if (add) {
      __builtin_prefetch(&books_[itch::load_be16(p + 1)]);
    }
    if (type == 'U' && len >= 27) {
      index_.prefetch(itch::load_be64(p + 19));
    }
    return;
  }
  if (add) {
    const Side_& v = books_[itch::load_be16(p + 1)].sides[p[19] == 'B' ? 0 : 1];
    if (!v.empty()) {
      __builtin_prefetch(&v.back());
      __builtin_prefetch(&levels_[v.back().level]);
    }
    return;
  }
  const std::uint32_t s = index_.find(ref);  // stage 0 brought the entry into cache
  if (s != kNil) {
    __builtin_prefetch(&slots_[s]);
  }
}

// ---- queries ----------------------------------------------------------------

std::optional<LevelInfo> FastBook::best(Locate locate, Side side) const {
  const Side_& v = books_[locate].sides[static_cast<int>(side)];
  if (v.empty()) {
    return std::nullopt;
  }
  const Level& b = levels_[v.back().level];
  return LevelInfo{b.price, b.shares, b.count};
}

std::vector<LevelInfo> FastBook::levels(Locate locate, Side side, std::size_t depth) const {
  const Side_& v = books_[locate].sides[static_cast<int>(side)];
  std::vector<LevelInfo> out;
  for (auto it = v.rbegin(); it != v.rend() && out.size() < depth; ++it) {
    const Level& lv = levels_[it->level];
    out.push_back({lv.price, lv.shares, lv.count});
  }
  return out;
}

std::optional<FastBook::Resting> FastBook::resting(OrderRef ref) const {
  const std::uint32_t s = find(ref);
  if (s == kNil) {
    return std::nullopt;
  }
  const Slot& o = slots_[s];
  return Resting{o.locate, o.side, o.price, o.shares};
}

bool FastBook::at_head(OrderRef ref) const {
  const std::uint32_t s = find(ref);
  return s != kNil && slots_[s].prev == kNil;
}

std::optional<FastBook::OrderInfo> FastBook::order(OrderRef ref) const {
  const std::uint32_t s = find(ref);
  if (s == kNil) {
    return std::nullopt;
  }
  const Slot& o = slots_[s];
  OrderInfo info{o.locate, o.side, o.price, o.shares, 0, 0};
  for (std::uint32_t p = o.prev; p != kNil; p = slots_[p].prev) {
    ++info.rank;
    info.shares_ahead += slots_[p].shares;
  }
  return info;
}

std::uint64_t FastBook::checksum() const {
  std::uint64_t sum = 0;
  for (std::size_t l = 0; l < books_.size(); ++l) {
    for (const Side side : {Side::buy, Side::sell}) {
      const Side_& v = books_[l].sides[static_cast<int>(side)];
      for (const PriceRef& pr : v) {  // the sum is order-independent across levels
        const Level& lv = levels_[pr.level];
        sum += mix({l, static_cast<std::uint64_t>(side), lv.price, lv.shares, lv.count});
        std::uint64_t rank = 0;
        for (std::uint32_t s = lv.head; s != kNil; s = slots_[s].next) {
          sum += mix({slots_[s].ref, l, static_cast<std::uint64_t>(side), lv.price, slots_[s].shares, rank++});
        }
      }
    }
  }
  return sum;
}

}  // namespace itchbook::book
