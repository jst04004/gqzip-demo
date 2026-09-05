#pragma once

#include "common.hpp"
#include <vector>
#include <string_view>
#include <span>

namespace gqzip {

struct RansFreqTable {
    static constexpr size_t SCALE_BITS = 12;
    static constexpr uint32_t TOTAL_FREQ = 1 << SCALE_BITS; // 4096

    std::array<uint16_t, 256> freqs{};
    std::array<uint16_t, 257> cum_freqs{};

    void build_from_data(const uint8_t* data, size_t size);
};

class RansCodec {
public:
    // Encodes a byte buffer into an rANS compressed bitstream
    static std::vector<uint8_t> encode(const uint8_t* data, size_t size);

    // Decodes an rANS bitstream back to the original byte buffer
    static std::vector<uint8_t> decode(const uint8_t* compressed, size_t compressed_size, size_t original_size);

    // High-performance block entropy compression for columnar streams
    static std::vector<uint8_t> compress_stream(const uint8_t* data, size_t size);
    static std::vector<uint8_t> decompress_stream(const uint8_t* compressed, size_t compressed_size, size_t original_size);

    // Order-1 Markov Context-Conditioned rANS Encoder/Decoder P(R_i | R_{i-1}, Q)
    static std::vector<uint8_t> encode_order1(const uint8_t* data, size_t size);
    static std::vector<uint8_t> decode_order1(const uint8_t* compressed, size_t compressed_size, size_t original_size);

    // Optimized Lossless Dual-Stream Residual Coder (Claim 3 ZigZag + Order-1 Zero-Run rANS)
    static std::vector<uint8_t> compress_residuals(std::span<const int8_t> residuals);
    static std::vector<int8_t> decompress_residuals(const uint8_t* compressed, size_t compressed_size, size_t original_size);
};

} // namespace gqzip
