#pragma once

#include "common.hpp"
#include "fastq_record.hpp"
#include "quantizer.hpp"
#include "rans_codec.hpp"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

namespace gqzip {

struct BlockData {
    uint32_t record_count = 0;
    uint32_t uncompressed_bytes = 0;
    uint16_t flags = 0; // Bit 0: has_residuals
    BinningLevel binning_level = BinningLevel::LEVEL_3_ADAPTIVE_CONTEXT;

    std::vector<uint8_t> header_stream;
    std::vector<uint8_t> dna_stream;
    std::vector<uint8_t> quality_stream;
    std::vector<uint8_t> residual_stream;

    uint32_t crc32 = 0;

    void serialize(std::ostream& out) const {
        // Write block magic / markers
        uint32_t marker = 0x4B4C4247; // 'GBLK'
        out.write(reinterpret_cast<const char*>(&marker), 4);
        out.write(reinterpret_cast<const char*>(&record_count), 4);
        out.write(reinterpret_cast<const char*>(&uncompressed_bytes), 4);
        out.write(reinterpret_cast<const char*>(&flags), 2);
        uint8_t b_level = static_cast<uint8_t>(binning_level);
        out.write(reinterpret_cast<const char*>(&b_level), 1);

        // Stream sizes
        uint32_t h_size = static_cast<uint32_t>(header_stream.size());
        uint32_t d_size = static_cast<uint32_t>(dna_stream.size());
        uint32_t q_size = static_cast<uint32_t>(quality_stream.size());
        uint32_t r_size = static_cast<uint32_t>(residual_stream.size());

        out.write(reinterpret_cast<const char*>(&h_size), 4);
        out.write(reinterpret_cast<const char*>(&d_size), 4);
        out.write(reinterpret_cast<const char*>(&q_size), 4);
        out.write(reinterpret_cast<const char*>(&r_size), 4);

        if (h_size > 0) out.write(reinterpret_cast<const char*>(header_stream.data()), h_size);
        if (d_size > 0) out.write(reinterpret_cast<const char*>(dna_stream.data()), d_size);
        if (q_size > 0) out.write(reinterpret_cast<const char*>(quality_stream.data()), q_size);
        if (r_size > 0) out.write(reinterpret_cast<const char*>(residual_stream.data()), r_size);

        out.write(reinterpret_cast<const char*>(&crc32), 4);
    }

    bool deserialize(std::istream& in) {
        uint32_t marker = 0;
        if (!in.read(reinterpret_cast<char*>(&marker), 4) || marker != 0x4B4C4247) {
            return false;
        }

        in.read(reinterpret_cast<char*>(&record_count), 4);
        in.read(reinterpret_cast<char*>(&uncompressed_bytes), 4);
        in.read(reinterpret_cast<char*>(&flags), 2);
        uint8_t b_level = 3;
        in.read(reinterpret_cast<char*>(&b_level), 1);
        binning_level = static_cast<BinningLevel>(b_level);

        uint32_t h_size = 0, d_size = 0, q_size = 0, r_size = 0;
        in.read(reinterpret_cast<char*>(&h_size), 4);
        in.read(reinterpret_cast<char*>(&d_size), 4);
        in.read(reinterpret_cast<char*>(&q_size), 4);
        in.read(reinterpret_cast<char*>(&r_size), 4);

        header_stream.resize(h_size);
        dna_stream.resize(d_size);
        quality_stream.resize(q_size);
        residual_stream.resize(r_size);

        if (h_size > 0) in.read(reinterpret_cast<char*>(header_stream.data()), h_size);
        if (d_size > 0) in.read(reinterpret_cast<char*>(dna_stream.data()), d_size);
        if (q_size > 0) in.read(reinterpret_cast<char*>(quality_stream.data()), q_size);
        if (r_size > 0) in.read(reinterpret_cast<char*>(residual_stream.data()), r_size);

        in.read(reinterpret_cast<char*>(&crc32), 4);
        return true;
    }
};

} // namespace gqzip
