#include "itchbook/book/fast_book.hpp"

#include <bit>
#include <limits>

namespace itchbook::book {

FastBook::FastBook(std::size_t expected_orders) : books_(std::size_t{std::numeric_limits<Locate>::max()} + 1) {
  // Keep the index at most half full: size = next power of two >= 2 * expected.
  const std::size_t size = std::bit_ceil(std::max<std::size_t>(expected_orders * 2, 1024));
  index_.assign(size, Bucket{0, kNil});
  mask_ = size - 1;
  shift_ = static_cast<unsigned>(64 - std::countr_zero(size));
  slots_.reserve(expected_orders);
}

// ---- order index ---------------------------------------------------------------

std::uint32_t FastBook::find(OrderRef ref) const noexcept {
  for (std::size_t i = home(ref);; i = (i + 1) & mask_) {
    const Bucket& b = index_[i];
    if (b.slot == kNil) {
      return kNil;
    }
    if (b.ref == ref) {
      return b.slot;
    }
  }
}

void FastBook::index_insert(OrderRef ref, std::uint32_t slot) {
  if ((live_ + 1) * 2 > index_.size()) {
    grow_index();
  }
  std::size_t i = home(ref);
  while (index_[i].slot != kNil) {
    i = (i + 1) & mask_;
  }
  index_[i] = {ref, slot};
}

void FastBook::index_erase(OrderRef ref) {
  std::size_t i = home(ref);
  while (index_[i].ref != ref || index_[i].slot == kNil) {
    i = (i + 1) & mask_;
  }
  // Backward-shift deletion: pull later entries of the same probe run into the
  // hole when that brings them no further from their home bucket than before.
  for (std::size_t j = (i + 1) & mask_;; j = (j + 1) & mask_) {
    if (index_[j].slot == kNil) {
      break;
    }
    const std::size_t dist_home = (j - home(index_[j].ref)) & mask_;
    const std::size_t dist_hole = (j - i) & mask_;
    if (dist_home >= dist_hole) {
      index_[i] = index_[j];
      i = j;
    }
  }
  index_[i].slot = kNil;
}

void FastBook::grow_index() {
  const std::vector<Bucket> old = std::move(index_);
  index_.assign(old.size() * 2, Bucket{0, kNil});
  mask_ = index_.size() - 1;
  --shift_;
  for (const Bucket& b : old) {
    if (b.slot != kNil) {
      std::size_t i = home(b.ref);
      while (index_[i].slot != kNil) {
        i = (i + 1) & mask_;
      }
      index_[i] = b;
    }
  }
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

FastBook::Level* FastBook::find_level(Locate l, Side s, Price p) noexcept {
  Side_& v = books_[l].sides[static_cast<int>(s)];
  const std::size_t i = lower(v, p, [s](Price a, Price b) { return better(s, a, b); });
  return (i < v.size() && v[i].price == p) ? &v[i] : nullptr;
}

const FastBook::Level* FastBook::find_level(Locate l, Side s, Price p) const noexcept {
  const Side_& v = books_[l].sides[static_cast<int>(s)];
  const std::size_t i = lower(v, p, [s](Price a, Price b) { return better(s, a, b); });
  return (i < v.size() && v[i].price == p) ? &v[i] : nullptr;
}

FastBook::Level& FastBook::level_for_insert(Locate l, Side s, Price p) {
  Side_& v = books_[l].sides[static_cast<int>(s)];
  const std::size_t i = lower(v, p, [s](Price a, Price b) { return better(s, a, b); });
  if (i < v.size() && v[i].price == p) {
    return v[i];
  }
  ++level_count_;
  return *v.insert(v.begin() + static_cast<std::ptrdiff_t>(i), Level{p, 0, 0, kNil, kNil});
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
  Level& lvl = level_for_insert(locate, side, price);
  slots_[s] = Slot{ref, price, shares, lvl.tail, kNil, locate, side};
  if (lvl.tail != kNil) {
    slots_[lvl.tail].next = s;
  } else {
    lvl.head = s;
  }
  lvl.tail = s;
  ++lvl.count;
  lvl.shares += shares;
  index_insert(ref, s);
  ++live_;
  stats_.max_live_orders = std::max<std::uint64_t>(stats_.max_live_orders, live_);
}

void FastBook::erase(std::uint32_t s) {
  const Slot o = slots_[s];
  Side_& v = books_[o.locate].sides[static_cast<int>(o.side)];
  Level* lvl = find_level(o.locate, o.side, o.price);
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
    v.erase(v.begin() + (lvl - v.data()));
    --level_count_;
  }
  index_erase(o.ref);
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
  find_level(o.locate, o.side, o.price)->shares -= n;
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

// ---- queries ----------------------------------------------------------------

std::optional<LevelInfo> FastBook::best(Locate locate, Side side) const {
  const Side_& v = books_[locate].sides[static_cast<int>(side)];
  if (v.empty()) {
    return std::nullopt;
  }
  return LevelInfo{v.back().price, v.back().shares, v.back().count};
}

std::vector<LevelInfo> FastBook::levels(Locate locate, Side side, std::size_t depth) const {
  const Side_& v = books_[locate].sides[static_cast<int>(side)];
  std::vector<LevelInfo> out;
  for (auto it = v.rbegin(); it != v.rend() && out.size() < depth; ++it) {
    out.push_back({it->price, it->shares, it->count});
  }
  return out;
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
      for (const Level& lv : v) {  // the sum is order-independent across levels
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
