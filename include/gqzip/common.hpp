#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <span>
#include <stdexcept>
#include <array>
#include <algorithm>
#include <cstring>

namespace gqzip {

// Magic identification bytes: "GQZ\x01" in little-endian (0x015A5147)
inline constexpr uint32_t MAGIC_HEADER = 0x015A5147;
inline constexpr uint16_t FORMAT_VERSION = 1;
inline constexpr size_t DEFAULT_BLOCK_RECORDS = 50000;
inline constexpr size_t MAX_READ_LENGTH = 16384;
inline constexpr size_t BOUNDED_RING_BUFFER_SLOTS = 8;
inline constexpr size_t MAX_RING_BUFFER_BYTES = 128 * 1024 * 1024; // 128 MB

enum class BinningLevel : uint8_t {
    LEVEL_1_ILLUMINA8 = 1,      // Standard Illumina 8-bin scheme
    LEVEL_2_COARSE4 = 2,        // Coarse 4-bin scheme
    LEVEL_3_ADAPTIVE_CONTEXT = 3,// Context-aware 1D window dynamic quantizer (Default)
    LEVEL_4_BINARY = 4,         // Binary high/low quality thresholding
    LEVEL_5_LOSSLESS = 5        // 100% Bit-exact lossless restoration with residual stream
};

enum class CompressionMode : uint8_t {
    COMPRESS = 0,
    DECOMPRESS = 1
};

struct CompressionOptions {
    BinningLevel binning_level = BinningLevel::LEVEL_3_ADAPTIVE_CONTEXT;
    bool lossless = false;
    size_t num_threads = 4;
    size_t block_records = DEFAULT_BLOCK_RECORDS;
    bool verbose = false;
};

inline constexpr std::array<uint32_t, 256> generate_crc32_table() {
    std::array<uint32_t, 256> table{};
    for (uint32_t i = 0; i < 256; ++i) {
        uint32_t crc = i;
        for (uint32_t j = 0; j < 8; ++j) {
            crc = (crc & 1) ? ((crc >> 1) ^ 0xEDB88320U) : (crc >> 1);
        }
        table[i] = crc;
    }
    return table;
}
inline constexpr std::array<uint32_t, 256> CRC32_TABLE = generate_crc32_table();

// Fast CRC32-C implementation table
class CRC32 {
public:
    static uint32_t calculate(const uint8_t* data, size_t length, uint32_t initial = 0xFFFFFFFF) {
        uint32_t crc = initial;
        for (size_t i = 0; i < length; ++i) {
            uint8_t index = static_cast<uint8_t>((crc ^ data[i]) & 0xFF);
            crc = (crc >> 8) ^ CRC32_TABLE[index];
        }
        return crc ^ 0xFFFFFFFF;
    }
};

} // namespace gqzip
