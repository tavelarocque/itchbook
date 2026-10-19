#pragma once
// Byte sources with bounded memory. A full NASDAQ day is 4-18 GB gzipped and
// 10-40 GB raw, far more than the 8 GB dev laptop has, so nothing here ever
// holds more than a fixed-size buffer.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <zlib.h>

namespace itchbook::io {

class ByteSource {
 public:
  ByteSource() = default;
  ByteSource(const ByteSource&) = delete;
  ByteSource& operator=(const ByteSource&) = delete;
  ByteSource(ByteSource&&) = delete;
  ByteSource& operator=(ByteSource&&) = delete;
  virtual ~ByteSource() = default;

  // Copy up to `cap` bytes into `dst`. Returns 0 only at end of input; may
  // return fewer than `cap` bytes at any time. Throws on I/O or format errors.
  virtual std::size_t read(std::uint8_t* dst, std::size_t cap) = 0;

  // True if the input ended in the middle of a compressed stream (for example
  // a partial download). Only meaningful after read() has returned 0.
  [[nodiscard]] virtual bool truncated() const noexcept { return false; }
};

// Plain file, or stdin when path is "-".
class FileSource final : public ByteSource {
 public:
  explicit FileSource(const std::string& path);
  ~FileSource() override;
  std::size_t read(std::uint8_t* dst, std::size_t cap) override;

  // Read up to n bytes ahead without consuming them (used to sniff gzip magic).
  std::span<const std::uint8_t> peek(std::size_t n);

 private:
  std::FILE* f_ = nullptr;
  bool owned_ = false;
  std::vector<std::uint8_t> pending_;
  std::size_t pending_pos_ = 0;
};

// Streaming gzip decoder. Handles multi-member files (concatenated gzip
// streams) and reports a truncated final member instead of throwing, so a
// partial download can still be used for development.
class GzipSource final : public ByteSource {
 public:
  explicit GzipSource(std::unique_ptr<ByteSource> inner, std::size_t in_buffer_bytes = 1U << 18U);
  ~GzipSource() override;
  std::size_t read(std::uint8_t* dst, std::size_t cap) override;
  [[nodiscard]] bool truncated() const noexcept override { return truncated_; }
  [[nodiscard]] std::uint64_t compressed_bytes_read() const noexcept { return compressed_read_; }

 private:
  std::unique_ptr<ByteSource> inner_;
  std::vector<std::uint8_t> in_;
  z_stream zs_{};
  bool in_eof_ = false;
  bool in_member_ = false;  // consumed header bytes of a member not yet finished
  bool done_ = false;
  bool truncated_ = false;
  std::uint64_t compressed_read_ = 0;
};

// In-memory bytes, handed out at most `max_chunk` at a time. Tests use small
// chunks to force frames across every possible buffer boundary.
class MemorySource final : public ByteSource {
 public:
  explicit MemorySource(std::span<const std::uint8_t> bytes, std::size_t max_chunk = SIZE_MAX)
      : bytes_(bytes), max_chunk_(max_chunk == 0 ? 1 : max_chunk) {}
  std::size_t read(std::uint8_t* dst, std::size_t cap) override;

 private:
  std::span<const std::uint8_t> bytes_;
  std::size_t pos_ = 0;
  std::size_t max_chunk_;
};

// Open a file (or "-" for stdin), decompressing transparently if it starts
// with the gzip magic bytes 1f 8b.
std::unique_ptr<ByteSource> open_source(const std::string& path);

}  // namespace itchbook::io
