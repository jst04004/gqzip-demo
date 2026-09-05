#pragma once

#include "common.hpp"
#include "fastq_record.hpp"
#include "entropy_model.hpp"
#include <string>
#include <vector>
#include <span>
#include <cstdint>

namespace gqzip {

struct QuantizedResult {
    std::string quantized_quality;
    alignas(32) std::vector<int8_t> residuals; // Empty if not in lossless mode (Claim 3 signed diff stream)
};

class Quantizer {
public:
    explicit Quantizer(BinningLevel level = BinningLevel::LEVEL_3_ADAPTIVE_CONTEXT, bool lossless = false)
        : m_level(level), m_lossless(lossless) {}

    // Bitwise ZigZag Mapping (Claim 3): Z(x) = (x << 1) ^ (x >> 31)
    static inline uint8_t zigzag_encode(int8_t diff) noexcept {
        int32_t val = static_cast<int32_t>(diff);
        return static_cast<uint8_t>((val << 1) ^ (val >> 31));
    }

    static inline int8_t zigzag_decode(uint8_t z) noexcept {
        uint32_t u = static_cast<uint32_t>(z);
        return static_cast<int8_t>((u >> 1) ^ -(u & 1));
    }

    // Quantizes a single Phred quality string given its corresponding sequence
    QuantizedResult quantize(std::string_view sequence, std::string_view quality) const;

    // Quantizes an entire columnar batch
    std::vector<QuantizedResult> quantize_batch(const FastqBatch& batch) const;

    // Restores original Phred quality string from quantized string and residual stream (C++20 std::span)
    static std::string restore(std::string_view quantized_quality, std::span<const int8_t> residuals);

    // Static binning lookup methods
    static uint8_t quantize_illumina8(uint8_t q);
    static uint8_t quantize_coarse4(uint8_t q);
    static uint8_t quantize_binary2(uint8_t q);
    static uint8_t quantize_adaptive_context(uint8_t q, float local_entropy, float cycle_factor);

private:
    BinningLevel m_level;
    bool m_lossless;
};

} // namespace gqzip

