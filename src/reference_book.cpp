#include "itchbook/book/reference_book.hpp"

#include <algorithm>
#include <limits>

namespace itchbook::book {

ReferenceBook::ReferenceBook() : books_(std::size_t{std::numeric_limits<Locate>::max()} + 1) {
  orders_.reserve(std::size_t{1} << 22U);
}

void ReferenceBook::check_locate(Locate header, const Order& o) noexcept {
  if (header != o.locate) {
    ++stats_.locate_mismatch;
  }
}

void ReferenceBook::insert(OrderRef ref, Locate locate, Side side, Shares shares, Price price) {
  SymbolBook& b = books_[locate];
  with_side(b, side, [&](auto& m) {
    auto [lit, created] = m.try_emplace(price);
    if (created) {
      ++level_count_;
    }
    Level& lvl = lit->second;
    lvl.shares += shares;
    lvl.queue.push_back(ref);
    orders_.emplace(ref, Order{locate, side, price, shares, std::prev(lvl.queue.end())});
  });
  stats_.max_live_orders = std::max<std::uint64_t>(stats_.max_live_orders, orders_.size());
}

void ReferenceBook::erase(std::unordered_map<OrderRef, Order>::iterator it) {
  const Order& o = it->second;
  SymbolBook& b = books_[o.locate];
  with_side(b, o.side, [&](auto& m) {
    auto lit = m.find(o.price);
    Level& lvl = lit->second;
    lvl.shares -= o.shares;
    lvl.queue.erase(o.pos);
    if (lvl.queue.empty()) {
      m.erase(lit);
      --level_count_;
    }
  });
  orders_.erase(it);
}

void ReferenceBook::reduce(std::unordered_map<OrderRef, Order>::iterator it, Shares n) {
  Order& o = it->second;
  if (n >= o.shares) {
    erase(it);
    return;
  }
  with_side(books_[o.locate], o.side, [&](auto& m) { m.find(o.price)->second.shares -= n; });
  o.shares -= n;
}

void ReferenceBook::add(Locate locate, OrderRef ref, char side, Shares shares, Price price) {
  if (side != 'B' && side != 'S') {
    ++stats_.bad_side;
    return;
  }
  if (orders_.contains(ref)) {
    ++stats_.duplicate_add;
    return;
  }
  ++stats_.adds;
  insert(ref, locate, side == 'B' ? Side::buy : Side::sell, shares, price);
}

void ReferenceBook::execute(Locate locate, OrderRef ref, Shares shares) {
  auto it = orders_.find(ref);
  if (it == orders_.end()) {
    ++stats_.orphan_execute;
    return;
  }
  ++stats_.executes;
  const Order& o = it->second;
  check_locate(locate, o);
  const bool at_head =
      with_side(books_[o.locate], o.side, [&](auto& m) { return m.find(o.price)->second.queue.begin() == o.pos; });
  ++(at_head ? stats_.executed_at_head : stats_.executed_behind_head);
  if (shares > o.shares) {
    ++stats_.over_execute;
  }
  reduce(it, shares);
}

void ReferenceBook::cancel(Locate locate, OrderRef ref, Shares shares) {
  auto it = orders_.find(ref);
  if (it == orders_.end()) {
    ++stats_.orphan_cancel;
    return;
  }
  ++stats_.cancels;
  check_locate(locate, it->second);
  if (shares > it->second.shares) {
    ++stats_.over_cancel;
  }
  if (shares >= it->second.shares) {
    ++stats_.cancel_to_zero;
  }
  reduce(it, shares);
}

void ReferenceBook::remove(Locate locate, OrderRef ref) {
  auto it = orders_.find(ref);
  if (it == orders_.end()) {
    ++stats_.orphan_delete;
    return;
  }
  ++stats_.deletes;
  check_locate(locate, it->second);
  erase(it);
}

void ReferenceBook::replace(Locate locate, OrderRef old_ref, OrderRef new_ref, Shares shares, Price price) {
  auto it = orders_.find(old_ref);
  if (it == orders_.end()) {
    ++stats_.orphan_replace;
    return;
  }
  if (new_ref != old_ref && orders_.contains(new_ref)) {
    ++stats_.duplicate_replace;
    return;
  }
  ++stats_.replaces;
  check_locate(locate, it->second);
  const Locate l = it->second.locate;
  const Side s = it->second.side;
  erase(it);
  if (shares == 0) {
    ++stats_.replace_to_zero;
    return;
  }
  insert(new_ref, l, s, shares, price);
}

std::optional<LevelInfo> ReferenceBook::best(Locate locate, Side side) const {
  return with_side(books_[locate], side, [](const auto& m) -> std::optional<LevelInfo> {
    if (m.empty()) {
      return std::nullopt;
    }
    const auto& [p, lvl] = *m.begin();
    return LevelInfo{p, lvl.shares, static_cast<std::uint32_t>(lvl.queue.size())};
  });
}

std::vector<LevelInfo> ReferenceBook::levels(Locate locate, Side side, std::size_t depth) const {
  std::vector<LevelInfo> out;
  with_side(books_[locate], side, [&](const auto& m) {
    for (const auto& [p, lvl] : m) {
      if (out.size() == depth) {
        break;
      }
      out.push_back({p, lvl.shares, static_cast<std::uint32_t>(lvl.queue.size())});
    }
  });
  return out;
}

bool ReferenceBook::at_head(OrderRef ref) const {
  const auto it = orders_.find(ref);
  if (it == orders_.end()) {
    return false;
  }
  const Order& o = it->second;
  return with_side(books_[o.locate], o.side,
                   [&](const auto& m) { return m.find(o.price)->second.queue.begin() == o.pos; });
}

std::optional<ReferenceBook::OrderInfo> ReferenceBook::order(OrderRef ref) const {
  const auto it = orders_.find(ref);
  if (it == orders_.end()) {
    return std::nullopt;
  }
  const Order& o = it->second;
  OrderInfo info{o.locate, o.side, o.price, o.shares, 0, 0};
  with_side(books_[o.locate], o.side, [&](const auto& m) {
    for (const OrderRef r : m.find(o.price)->second.queue) {
      if (r == ref) {
        break;
      }
      ++info.rank;
      info.shares_ahead += orders_.at(r).shares;
    }
  });
  return info;
}

std::uint64_t ReferenceBook::checksum() const {
  std::uint64_t sum = 0;
  for (std::size_t l = 0; l < books_.size(); ++l) {
    for (const Side side : {Side::buy, Side::sell}) {
      with_side(books_[l], side, [&](const auto& m) {
        for (const auto& [p, lvl] : m) {
          sum += mix({l, static_cast<std::uint64_t>(side), p, lvl.shares, lvl.queue.size()});
          std::uint64_t rank = 0;
          for (const OrderRef r : lvl.queue) {
            sum += mix({r, l, static_cast<std::uint64_t>(side), p, orders_.at(r).shares, rank++});
          }
        }
      });
    }
  }
  return sum;
}

}  // namespace itchbook::book
