#pragma once

#include "common.hpp"
#include "fastq_record.hpp"
#include <string_view>
#include <vector>
#include <istream>

namespace gqzip {

class SIMDScanner {
public:
    // Scans a byte buffer for newline ('\n') positions using SIMD vector instructions
    // with portable SWAR (SIMD within a register) fallback.
    static std::vector<size_t> find_newlines(const char* data, size_t length);

    // Parses raw FASTQ text chunk into columnar FastqBatch with comprehensive validation
    static bool parse_fastq_chunk(
        std::string_view chunk,
        FastqBatch& batch,
        std::string& remaining_tail,
        bool is_eof = false
    );

    // Stream reader that reads from std::istream in bounded chunks
    static bool read_batch(
        std::istream& stream,
        FastqBatch& batch,
        size_t max_records = DEFAULT_BLOCK_RECORDS
    );
};

} // namespace gqzip
