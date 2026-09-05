#pragma once

#include "common.hpp"
#include "fastq_record.hpp"
#include "quantizer.hpp"
#include "simd_scanner.hpp"
#include "block_stream.hpp"
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>

namespace gqzip {

struct CompressionStats {
    uint64_t total_records = 0;
    uint64_t raw_fastq_bytes = 0;
    uint64_t compressed_bytes = 0;
    uint64_t raw_quality_bytes = 0;
    uint64_t compressed_quality_bytes = 0;
    double compression_time_sec = 0.0;
    double decompression_time_sec = 0.0;
    size_t peak_memory_bytes = 0;

    double compression_ratio() const {
        return (raw_fastq_bytes > 0) ? static_cast<double>(raw_fastq_bytes) / compressed_bytes : 1.0;
    }

    double quality_space_saving_pct() const {
        return (raw_quality_bytes > 0) ? (1.0 - static_cast<double>(compressed_quality_bytes) / raw_quality_bytes) * 100.0 : 0.0;
    }

    double throughput_mb_s() const {
        return (compression_time_sec > 0) ? (static_cast<double>(raw_fastq_bytes) / (1024.0 * 1024.0)) / compression_time_sec : 0.0;
    }
};

class CompressionEngine {
public:
    explicit CompressionEngine(const CompressionOptions& options = {})
        : m_options(options), m_quantizer(options.binning_level, options.lossless) {}

    // Compress standard FASTQ stream into binary .gqz format
    CompressionStats compress(std::istream& in, std::ostream& out);

    // Decompress binary .gqz stream into standard FASTQ format
    CompressionStats decompress(std::istream& in, std::ostream& out);

    // Block level compression/decompression helpers
    BlockData compress_batch(const FastqBatch& batch);
    FastqBatch decompress_block(const BlockData& block);

private:
    CompressionOptions m_options;
    Quantizer m_quantizer;

    // DNA 2-bit packing (A=0, C=1, G=2, T=3, with 'N' mask)
    static std::vector<uint8_t> pack_dna(const std::vector<std::string>& sequences);
    static std::vector<std::string> unpack_dna(const uint8_t* data, size_t size, size_t record_count);

    // Header tokenization & delta coding
    static std::vector<uint8_t> encode_headers(const std::vector<std::string>& headers);
    static std::vector<std::string> decode_headers(const uint8_t* data, size_t size, size_t record_count);
};

} // namespace gqzip
