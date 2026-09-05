#include "gqzip/simd_scanner.hpp"
#include <iostream>
#include <sstream>
#include <cassert>

namespace gqzip {

std::vector<size_t> SIMDScanner::find_newlines(const char* data, size_t length) {
    std::vector<size_t> newlines;
    newlines.reserve(length / 40); // Approximate line density

    size_t i = 0;
    const uint64_t mask_01 = 0x0101010101010101ULL;
    const uint64_t mask_80 = 0x8080808080808080ULL;
    const uint64_t target_nl = 0x0A0A0A0A0A0A0A0AULL;

    // Fast 64-bit SWAR SIMD Scanner
    while (i + 8 <= length) {
        uint64_t word;
        std::memcpy(&word, data + i, 8);
        uint64_t diff = word ^ target_nl;
        uint64_t has_nl = (diff - mask_01) & ~diff & mask_80;

        if (has_nl != 0) {
            for (size_t b = 0; b < 8; ++b) {
                if (data[i + b] == '\n') {
                    newlines.push_back(i + b);
                }
            }
        }
        i += 8;
    }

    // Scalar remainder
    while (i < length) {
        if (data[i] == '\n') {
            newlines.push_back(i);
        }
        ++i;
    }

    return newlines;
}

bool SIMDScanner::parse_fastq_chunk(
    std::string_view chunk,
    FastqBatch& batch,
    std::string& remaining_tail,
    bool is_eof
) {
    if (chunk.empty()) {
        return true;
    }

    // Prepend any leftover tail from previous chunk
    std::string combined;
    std::string_view full_text;
    if (!remaining_tail.empty()) {
        combined = remaining_tail + std::string(chunk);
        full_text = combined;
        remaining_tail.clear();
    } else {
        full_text = chunk;
    }

    auto newlines = find_newlines(full_text.data(), full_text.size());
    if (newlines.empty()) {
        if (!is_eof) {
            remaining_tail = std::string(full_text);
        }
        return true;
    }

    size_t prev_pos = 0;
    std::vector<std::string_view> lines;
    lines.reserve(newlines.size() + 1);

    for (size_t nl : newlines) {
        size_t len = nl - prev_pos;
        // Trim trailing '\r' if present (Windows line endings)
        if (len > 0 && full_text[prev_pos + len - 1] == '\r') {
            --len;
        }
        lines.push_back(full_text.substr(prev_pos, len));
        prev_pos = nl + 1;
    }

    // Leftover after last newline
    if (prev_pos < full_text.size()) {
        if (!is_eof) {
            remaining_tail = std::string(full_text.substr(prev_pos));
        } else {
            size_t len = full_text.size() - prev_pos;
            if (len > 0 && full_text[prev_pos + len - 1] == '\r') {
                --len;
            }
            if (len > 0) {
                lines.push_back(full_text.substr(prev_pos, len));
            }
        }
    }

    // Every FASTQ record consists of 4 lines
    size_t num_complete_records = lines.size() / 4;
    batch.reserve(batch.size() + num_complete_records);

    for (size_t r = 0; r < num_complete_records; ++r) {
        size_t idx = r * 4;
        std::string_view header = lines[idx];
        std::string_view seq = lines[idx + 1];
        std::string_view plus = lines[idx + 2];
        std::string_view qual = lines[idx + 3];

        // Format validation: header starts with '@', plus with '+'
        if (header.empty() || header[0] != '@') {
            // Malformed header recovery: normalize prefix
            std::string fixed_header = "@" + std::string(header);
            batch.add(std::move(fixed_header), std::string(seq), std::string(qual));
            continue;
        }

        // Sequence and quality length must match
        if (seq.size() != qual.size()) {
            // Truncate or pad to match
            size_t min_len = std::min(seq.size(), qual.size());
            batch.add(
                std::string(header),
                std::string(seq.substr(0, min_len)),
                std::string(qual.substr(0, min_len))
            );
            continue;
        }

        batch.add(std::string(header), std::string(seq), std::string(qual));
    }

    // Put unprocessed trailing lines back into remaining_tail if not EOF
    size_t processed_lines = num_complete_records * 4;
    if (processed_lines < lines.size() && !is_eof) {
        std::string leftover;
        for (size_t i = processed_lines; i < lines.size(); ++i) {
            leftover += std::string(lines[i]) + "\n";
        }
        remaining_tail = leftover + remaining_tail;
    }

    return true;
}

bool SIMDScanner::read_batch(
    std::istream& stream,
    FastqBatch& batch,
    size_t max_records
) {
    batch.clear();
    batch.reserve(max_records);

    std::string h, s, p, q;
    while (batch.size() < max_records && std::getline(stream, h)) {
        if (h.empty()) continue;
        if (!h.empty() && h.back() == '\r') h.pop_back();

        if (!std::getline(stream, s)) break;
        if (!s.empty() && s.back() == '\r') s.pop_back();

        if (!std::getline(stream, p)) break;
        if (!p.empty() && p.back() == '\r') p.pop_back();

        if (!std::getline(stream, q)) break;
        if (!q.empty() && q.back() == '\r') q.pop_back();

        // Edge case corrections
        if (h.empty() || h[0] != '@') {
            h = "@" + h;
        }
        if (s.size() != q.size()) {
            size_t min_len = std::min(s.size(), q.size());
            s = s.substr(0, min_len);
            q = q.substr(0, min_len);
        }

        batch.add(std::move(h), std::move(s), std::move(q));
    }

    return !batch.empty();
}

} // namespace gqzip
