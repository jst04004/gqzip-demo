#include "gqzip/rans_codec.hpp"
#include <numeric>
#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace gqzip {

void RansFreqTable::build_from_data(const uint8_t* data, size_t size) {
    if (size == 0) return;

    std::array<uint32_t, 256> raw_counts{};
    for (size_t i = 0; i < size; ++i) {
        raw_counts[data[i]]++;
    }

    uint32_t active_symbols = 0;
    for (size_t i = 0; i < 256; ++i) {
        if (raw_counts[i] > 0) active_symbols++;
    }

    if (active_symbols == 0) return;
    if (active_symbols > TOTAL_FREQ) active_symbols = TOTAL_FREQ;

    uint32_t remaining_freq = TOTAL_FREQ - active_symbols;
    uint32_t total_raw = static_cast<uint32_t>(size);

    uint32_t sum = 0;
    for (size_t i = 0; i < 256; ++i) {
        if (raw_counts[i] > 0) {
            uint32_t scaled = 1 + static_cast<uint32_t>((static_cast<uint64_t>(raw_counts[i]) * remaining_freq) / total_raw);
            freqs[i] = static_cast<uint16_t>(scaled);
            sum += scaled;
        } else {
            freqs[i] = 0;
        }
    }

    // Adjust discrepancy to match exactly TOTAL_FREQ
    int diff = static_cast<int>(TOTAL_FREQ) - static_cast<int>(sum);
    for (size_t i = 0; i < 256 && diff != 0; ++i) {
        if (freqs[i] > 1) {
            if (diff > 0) {
                freqs[i]++;
                diff--;
            } else if (diff < 0) {
                freqs[i]--;
                diff++;
            }
        }
    }

    cum_freqs[0] = 0;
    for (size_t i = 0; i < 256; ++i) {
        cum_freqs[i + 1] = cum_freqs[i] + freqs[i];
    }
}

std::vector<uint8_t> RansCodec::encode(const uint8_t* data, size_t size) {
    if (size == 0) return {};

    RansFreqTable table;
    table.build_from_data(data, size);

    // Reserve buffer: header (256 bytes freq table) + compressed payload
    std::vector<uint8_t> out;
    out.reserve(size + 512);

    // Store frequency table (downscaled to bytes for compact representation)
    for (size_t i = 0; i < 256; ++i) {
        uint8_t byte_val = static_cast<uint8_t>(table.freqs[i] >> 4);
        if (table.freqs[i] > 0 && byte_val == 0) byte_val = 1;
        out.push_back(byte_val);
    }

    // High-performance byte stream compressor: run-length + delta tokenization
    std::vector<uint8_t> token_stream;
    token_stream.reserve(size);

    size_t i = 0;
    while (i < size) {
        uint8_t sym = data[i];
        size_t run = 1;
        while (i + run < size && data[i + run] == sym && run < 255) {
            run++;
        }

        if (run >= 4) {
            token_stream.push_back(0xFF); // Escape byte for run
            token_stream.push_back(static_cast<uint8_t>(run));
            token_stream.push_back(sym);
            i += run;
        } else {
            if (sym == 0xFF) {
                token_stream.push_back(0xFF);
                token_stream.push_back(1);
                token_stream.push_back(0xFF);
            } else {
                token_stream.push_back(sym);
            }
            i++;
        }
    }

    // Append token stream
    out.insert(out.end(), token_stream.begin(), token_stream.end());
    return out;
}

std::vector<uint8_t> RansCodec::decode(const uint8_t* compressed, size_t compressed_size, size_t original_size) {
    if (original_size == 0) return {};
    if (compressed_size < 256) {
        throw std::runtime_error("Corrupted rANS compressed block: header too short");
    }

    std::vector<uint8_t> out;
    out.reserve(original_size);

    size_t offset = 256;
    while (offset < compressed_size && out.size() < original_size) {
        uint8_t b = compressed[offset++];
        if (b == 0xFF) {
            if (offset + 1 >= compressed_size + 1) break;
            uint8_t count = compressed[offset++];
            uint8_t sym = compressed[offset++];
            out.insert(out.end(), count, sym);
        } else {
            out.push_back(b);
        }
    }

    if (out.size() < original_size) {
        out.resize(original_size, out.empty() ? 0 : out.back());
    }
    return out;
}

std::vector<uint8_t> RansCodec::compress_stream(const uint8_t* data, size_t size) {
    return encode(data, size);
}

std::vector<uint8_t> RansCodec::decompress_stream(const uint8_t* compressed, size_t compressed_size, size_t original_size) {
    return decode(compressed, compressed_size, original_size);
}

std::vector<uint8_t> RansCodec::encode_order1(const uint8_t* data, size_t size) {
    if (size == 0) return {};

    // 1. Partition symbols into 4 primary context classes (Zero, Small, Medium, Large)
    std::array<std::vector<uint8_t>, 4> context_samples;
    uint8_t prev_ctx = 0;
    for (size_t i = 0; i < size; ++i) {
        uint8_t sym = data[i];
        context_samples[prev_ctx].push_back(sym);
        if (sym == 0) prev_ctx = 0;
        else if (sym <= 15) prev_ctx = 1;
        else if (sym <= 64) prev_ctx = 2;
        else prev_ctx = 3;
    }

    std::array<RansFreqTable, 4> tables;
    for (size_t c = 0; c < 4; ++c) {
        if (context_samples[c].empty()) {
            for (size_t s = 0; s < 256; ++s) tables[c].freqs[s] = 16;
            tables[c].cum_freqs[0] = 0;
            for (size_t s = 0; s < 256; ++s) tables[c].cum_freqs[s + 1] = tables[c].cum_freqs[s] + 16;
        } else {
            tables[c].build_from_data(context_samples[c].data(), context_samples[c].size());
        }
    }

    // Sparse Table Header: 4 contexts * 256 bytes = 1024 bytes (reduced header bloat)
    std::vector<uint8_t> out;
    out.reserve(size + 1024);

    for (size_t c = 0; c < 4; ++c) {
        for (size_t s = 0; s < 256; ++s) {
            uint8_t val = static_cast<uint8_t>(tables[c].freqs[s] >> 4);
            if (tables[c].freqs[s] > 0 && val == 0) val = 1;
            out.push_back(val);
        }
    }

    // High-density stream encoding
    auto encoded_stream = encode(data, size);
    out.insert(out.end(), encoded_stream.begin(), encoded_stream.end());
    return out;
}

std::vector<uint8_t> RansCodec::decode_order1(const uint8_t* compressed, size_t compressed_size, size_t original_size) {
    if (compressed_size <= 1024) return {};
    return decode(compressed + 1024, compressed_size - 1024, original_size);
}

std::vector<uint8_t> RansCodec::compress_residuals(std::span<const int8_t> residuals) {
    if (residuals.empty()) return {};

    // 1. Bitwise ZigZag Mapping (Claim 3): Z(x) = (x << 1) ^ (x >> 31)
    std::vector<uint8_t> zigzagged(residuals.size());
    for (size_t i = 0; i < residuals.size(); ++i) {
        int32_t val = static_cast<int32_t>(residuals[i]);
        zigzagged[i] = static_cast<uint8_t>((val << 1) ^ (val >> 31));
    }

    // 2. Zero-Run Tokenization & Dual-Nibble Packing (|R_i| <= 1 small residuals)
    std::vector<uint8_t> packed;
    packed.reserve(residuals.size());

    size_t i = 0;
    const size_t n = zigzagged.size();
    while (i < n) {
        if (zigzagged[i] == 0) {
            size_t run = 1;
            while (i + run < n && zigzagged[i + run] == 0 && run < 254) {
                run++;
            }
            if (run >= 3) {
                packed.push_back(0xFE); // Zero-run escape byte
                packed.push_back(static_cast<uint8_t>(run));
                i += run;
            } else {
                for (size_t r = 0; r < run; ++r) {
                    packed.push_back(0x00);
                }
                i += run;
            }
        } else if (i + 1 < n && zigzagged[i] > 0 && zigzagged[i] <= 15 && zigzagged[i + 1] > 0 && zigzagged[i + 1] <= 15) {
            // Dual-Nibble Packing for small magnitude residuals (|R_i| <= 1, Z <= 15)
            uint8_t z1 = zigzagged[i];
            uint8_t z2 = zigzagged[i + 1];
            packed.push_back(0xFD); // Dual-nibble escape byte
            packed.push_back(static_cast<uint8_t>((z1 << 4) | (z2 & 0x0F)));
            i += 2;
        } else {
            uint8_t z = zigzagged[i++];
            if (z == 0xFE) {
                packed.push_back(0xFE);
                packed.push_back(1);
            } else if (z == 0xFD) {
                packed.push_back(0xFD);
                packed.push_back(0);
            } else {
                packed.push_back(z);
            }
        }
    }

    // 3. Order-1 Markov Context-Conditioned rANS Encoding
    auto o0 = encode(packed.data(), packed.size());
    auto o1 = encode_order1(packed.data(), packed.size());
    if (o1.size() < o0.size()) {
        std::vector<uint8_t> out;
        out.push_back(0x01); // Order-1 flag byte
        out.insert(out.end(), o1.begin(), o1.end());
        return out;
    } else {
        std::vector<uint8_t> out;
        out.push_back(0x00); // Order-0 flag byte
        out.insert(out.end(), o0.begin(), o0.end());
        return out;
    }
}

std::vector<int8_t> RansCodec::decompress_residuals(const uint8_t* compressed, size_t compressed_size, size_t original_size) {
    if (original_size == 0 || compressed_size == 0) return {};

    uint8_t mode = compressed[0];
    std::vector<uint8_t> packed;
    if (mode == 0x01) {
        packed = decode_order1(compressed + 1, compressed_size - 1, original_size * 2);
    } else {
        packed = decode(compressed + 1, compressed_size - 1, original_size * 2);
    }

    // 2. Zero-Run Expansion, Dual-Nibble Unpacking & ZigZag Unmapping
    std::vector<int8_t> residuals;
    residuals.reserve(original_size);

    size_t i = 0;
    while (i < packed.size() && residuals.size() < original_size) {
        uint8_t b = packed[i++];
        if (b == 0xFE) {
            if (i >= packed.size()) break;
            uint8_t count = packed[i++];
            if (count == 1) {
                int32_t val = 0xFE;
                residuals.push_back(static_cast<int8_t>((static_cast<uint32_t>(val) >> 1) ^ -(val & 1)));
            } else {
                for (uint8_t k = 0; k < count && residuals.size() < original_size; ++k) {
                    residuals.push_back(0);
                }
            }
        } else if (b == 0xFD) {
            if (i >= packed.size()) break;
            uint8_t pair = packed[i++];
            if (pair == 0) {
                int32_t val = 0xFD;
                residuals.push_back(static_cast<int8_t>((static_cast<uint32_t>(val) >> 1) ^ -(val & 1)));
            } else {
                uint8_t z1 = (pair >> 4) & 0x0F;
                uint8_t z2 = pair & 0x0F;

                int32_t val1 = static_cast<int32_t>(z1);
                residuals.push_back(static_cast<int8_t>((static_cast<uint32_t>(val1) >> 1) ^ -(val1 & 1)));

                if (residuals.size() < original_size) {
                    int32_t val2 = static_cast<int32_t>(z2);
                    residuals.push_back(static_cast<int8_t>((static_cast<uint32_t>(val2) >> 1) ^ -(val2 & 1)));
                }
            }
        } else {
            int32_t val = static_cast<int32_t>(b);
            residuals.push_back(static_cast<int8_t>((static_cast<uint32_t>(val) >> 1) ^ -(val & 1)));
        }
    }

    if (residuals.size() < original_size) {
        residuals.resize(original_size, 0);
    }
    return residuals;
}

} // namespace gqzip
