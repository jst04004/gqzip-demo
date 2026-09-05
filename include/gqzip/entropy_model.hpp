#pragma once

#include "common.hpp"
#include <string_view>
#include <vector>
#include <cmath>
#include <array>

namespace gqzip {

class EntropyModel {
public:
    // Calculates Shannon entropy of nucleotide 4-base composition in window [i-w, i+w]
    static float calculate_local_entropy(
        std::string_view seq,
        size_t pos,
        size_t window_half_width = 3
    ) {
        size_t len = seq.size();
        if (len == 0) return 0.0f;

        size_t start = (pos >= window_half_width) ? (pos - window_half_width) : 0;
        size_t end = std::min(pos + window_half_width + 1, len);
        size_t win_len = end - start;
        if (win_len == 0) return 0.0f;

        std::array<uint32_t, 5> counts{}; // A, C, G, T, Other
        for (size_t k = start; k < end; ++k) {
            char b = seq[k];
            switch (b) {
                case 'A': case 'a': counts[0]++; break;
                case 'C': case 'c': counts[1]++; break;
                case 'G': case 'g': counts[2]++; break;
                case 'T': case 't': counts[3]++; break;
                default:            counts[4]++; break;
            }
        }

        float entropy = 0.0f;
        float inv_len = 1.0f / static_cast<float>(win_len);
        for (size_t b = 0; b < 4; ++b) {
            if (counts[b] > 0) {
                float p = static_cast<float>(counts[b]) * inv_len;
                entropy -= p * std::log2(p);
            }
        }

        return entropy;
    }

    // Calculates cycle position weighting: D(i) = 1.0 - alpha * (i / L)^gamma
    static float calculate_cycle_factor(
        size_t pos,
        size_t total_len,
        float alpha = 0.35f,
        float gamma = 1.8f
    ) {
        if (total_len <= 1) return 1.0f;
        float ratio = static_cast<float>(pos) / static_cast<float>(total_len);
        return 1.0f - alpha * std::pow(ratio, gamma);
    }
};

} // namespace gqzip
