#include "gqzip/quantizer.hpp"
#include <algorithm>
#include <cassert>

namespace gqzip {

uint8_t Quantizer::quantize_illumina8(uint8_t q) {
    if (q <= 1) return 0;
    if (q <= 9) return 6;
    if (q <= 19) return 15;
    if (q <= 24) return 22;
    if (q <= 29) return 27;
    if (q <= 34) return 33;
    if (q <= 39) return 37;
    return 40;
}

uint8_t Quantizer::quantize_coarse4(uint8_t q) {
    if (q <= 9) return 6;
    if (q <= 24) return 18;
    if (q <= 34) return 30;
    return 38;
}

uint8_t Quantizer::quantize_binary2(uint8_t q) {
    return (q < 20) ? 10 : 35;
}

uint8_t Quantizer::quantize_adaptive_context(uint8_t q, float local_entropy, float cycle_factor) {
    // Low-complexity / repetitive region: preserve fine-grained fidelity
    if (local_entropy < 1.0f) {
        return quantize_illumina8(q);
    }

    // High-complexity / diverse sequence: coarse binning eliminates noise without hurting mapping
    if (local_entropy >= 1.5f) {
        // Late-cycle degradation: slightly favor lower noise bin
        if (cycle_factor < 0.75f && q >= 28 && q <= 34) {
            return 30;
        }
        return quantize_coarse4(q);
    }

    // Intermediate 6-bin context:
    if (q <= 2) return 0;
    if (q <= 12) return 8;
    if (q <= 22) return 18;
    if (q <= 29) return 26;
    if (q <= 36) return 33;
    return 39;
}

QuantizedResult Quantizer::quantize(std::string_view sequence, std::string_view quality) const {
    size_t len = std::min(sequence.size(), quality.size());
    QuantizedResult result;
    result.quantized_quality.resize(len);

    if (m_lossless || m_level == BinningLevel::LEVEL_5_LOSSLESS) {
        result.residuals.resize(len);
    }

    for (size_t i = 0; i < len; ++i) {
        uint8_t raw_char = static_cast<uint8_t>(quality[i]);
        uint8_t raw_q = (raw_char >= 33) ? (raw_char - 33) : 0;
        uint8_t quant_q = raw_q;

        switch (m_level) {
            case BinningLevel::LEVEL_1_ILLUMINA8:
                quant_q = quantize_illumina8(raw_q);
                break;
            case BinningLevel::LEVEL_2_COARSE4:
                quant_q = quantize_coarse4(raw_q);
                break;
            case BinningLevel::LEVEL_3_ADAPTIVE_CONTEXT: {
                float entropy = EntropyModel::calculate_local_entropy(sequence, i, 3);
                float cycle = EntropyModel::calculate_cycle_factor(i, len);
                quant_q = quantize_adaptive_context(raw_q, entropy, cycle);
                break;
            }
            case BinningLevel::LEVEL_4_BINARY:
                quant_q = quantize_binary2(raw_q);
                break;
            case BinningLevel::LEVEL_5_LOSSLESS: {
                float entropy = EntropyModel::calculate_local_entropy(sequence, i, 3);
                float cycle = EntropyModel::calculate_cycle_factor(i, len);
                quant_q = quantize_adaptive_context(raw_q, entropy, cycle);
                break;
            }
        }

        // Clamp quantized Q to valid Phred range [0, 93] -> ASCII [33, 126]
        quant_q = std::min<uint8_t>(quant_q, 93);
        result.quantized_quality[i] = static_cast<char>(quant_q + 33);

        if (m_lossless || m_level == BinningLevel::LEVEL_5_LOSSLESS) {
            int diff = static_cast<int>(raw_q) - static_cast<int>(quant_q);
            result.residuals[i] = static_cast<int8_t>(diff);
        }
    }

    return result;
}

std::vector<QuantizedResult> Quantizer::quantize_batch(const FastqBatch& batch) const {
    std::vector<QuantizedResult> results;
    results.reserve(batch.size());
    for (size_t i = 0; i < batch.size(); ++i) {
        results.push_back(quantize(batch.sequences[i], batch.qualities[i]));
    }
    return results;
}

std::string Quantizer::restore(std::string_view quantized_quality, std::span<const int8_t> residuals) {
    if (residuals.empty()) {
        return std::string(quantized_quality);
    }

    size_t len = std::min(quantized_quality.size(), residuals.size());
    std::string raw(len, '!');
    for (size_t i = 0; i < len; ++i) {
        uint8_t quant_char = static_cast<uint8_t>(quantized_quality[i]);
        int quant_q = (quant_char >= 33) ? (quant_char - 33) : 0;
        int restored_q = quant_q + static_cast<int>(residuals[i]);
        restored_q = std::max(0, std::min(93, restored_q));
        raw[i] = static_cast<char>(restored_q + 33);
    }
    return raw;
}

} // namespace gqzip
