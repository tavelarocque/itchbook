#include "itchbook/io/source.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace itchbook::io {

// ---- FileSource ------------------------------------------------------------

FileSource::FileSource(const std::string& path) {
  if (path == "-") {
    f_ = stdin;
    return;
  }
  f_ = std::fopen(path.c_str(), "rb");
  if (f_ == nullptr) {
    throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
  }
  owned_ = true;
  // Large stdio buffer: fewer syscalls on multi-GB sequential reads.
  std::setvbuf(f_, nullptr, _IOFBF, 1U << 20U);
}

FileSource::~FileSource() {
  if (owned_ && f_ != nullptr) {
    std::fclose(f_);
  }
}

std::span<const std::uint8_t> FileSource::peek(std::size_t n) {
  while (pending_.size() - pending_pos_ < n) {
    std::uint8_t tmp[64];
    const std::size_t want = std::min(sizeof tmp, n - (pending_.size() - pending_pos_));
    const std::size_t got = std::fread(tmp, 1, want, f_);
    if (got == 0) {
      break;
    }
    pending_.insert(pending_.end(), tmp, tmp + got);
  }
  const auto avail = std::min(n, pending_.size() - pending_pos_);
  return {pending_.data() + pending_pos_, avail};
}

std::size_t FileSource::read(std::uint8_t* dst, std::size_t cap) {
  std::size_t out = 0;
  if (pending_pos_ < pending_.size()) {
    out = std::min(cap, pending_.size() - pending_pos_);
    std::memcpy(dst, pending_.data() + pending_pos_, out);
    pending_pos_ += out;
    if (pending_pos_ == pending_.size()) {
      pending_.clear();
      pending_pos_ = 0;
    }
    if (out == cap) {
      return out;
    }
  }
  const std::size_t got = std::fread(dst + out, 1, cap - out, f_);
  if (got < cap - out && std::ferror(f_) != 0) {
    throw std::runtime_error(std::string("read error: ") + std::strerror(errno));
  }
  return out + got;
}

// ---- GzipSource ------------------------------------------------------------

GzipSource::GzipSource(std::unique_ptr<ByteSource> inner, std::size_t in_buffer_bytes)
    : inner_(std::move(inner)), in_(in_buffer_bytes) {
  // 15 window bits + 32: accept a gzip or zlib header, detected automatically.
  if (inflateInit2(&zs_, 15 + 32) != Z_OK) {
    throw std::runtime_error("inflateInit2 failed");
  }
}

GzipSource::~GzipSource() { inflateEnd(&zs_); }

std::size_t GzipSource::read(std::uint8_t* dst, std::size_t cap) {
  if (done_ || cap == 0) {
    return 0;
  }
  constexpr std::size_t kMaxAvail = std::numeric_limits<uInt>::max();
  zs_.next_out = dst;
  zs_.avail_out = static_cast<uInt>(std::min(cap, kMaxAvail));
  const uInt out_cap = zs_.avail_out;

  while (zs_.avail_out > 0) {
    if (zs_.avail_in == 0) {
      if (in_eof_) {
        break;
      }
      const std::size_t n = inner_->read(in_.data(), in_.size());
      if (n == 0) {
        in_eof_ = true;
        break;
      }
      compressed_read_ += n;
      zs_.next_in = in_.data();
      zs_.avail_in = static_cast<uInt>(n);
    }
    in_member_ = true;
    const int rc = inflate(&zs_, Z_NO_FLUSH);
    if (rc == Z_STREAM_END) {
      // End of one gzip member. Another may follow (multi-member gzip).
      in_member_ = false;
      if (inflateReset(&zs_) != Z_OK) {
        throw std::runtime_error("inflateReset failed");
      }
      continue;
    }
    if (rc == Z_OK || rc == Z_BUF_ERROR) {
      continue;  // Z_BUF_ERROR: needs more input or output space; loop handles both.
    }
    throw std::runtime_error(std::string("gzip decode error: ") + (zs_.msg != nullptr ? zs_.msg : "unknown"));
  }

  const std::size_t produced = out_cap - zs_.avail_out;
  if (produced == 0 && in_eof_ && zs_.avail_in == 0) {
    done_ = true;
    // A member was started but never reached Z_STREAM_END: partial file. A
    // zero-length tail after the last member is fine.
    truncated_ = in_member_;
  }
  return produced;
}

// ---- MemorySource ----------------------------------------------------------

std::size_t MemorySource::read(std::uint8_t* dst, std::size_t cap) {
  const std::size_t n = std::min({cap, max_chunk_, bytes_.size() - pos_});
  if (n > 0) {
    std::memcpy(dst, bytes_.data() + pos_, n);
    pos_ += n;
  }
  return n;
}

// ---- open_source -----------------------------------------------------------

std::unique_ptr<ByteSource> open_source(const std::string& path) {
  auto file = std::make_unique<FileSource>(path);
  const auto head = file->peek(2);
  if (head.size() == 2 && head[0] == 0x1f && head[1] == 0x8b) {
    return std::make_unique<GzipSource>(std::move(file));
  }
  return file;
}

}  // namespace itchbook::io
