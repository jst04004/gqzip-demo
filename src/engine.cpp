#include "gqzip/engine.hpp"
#include <chrono>
#include <sstream>

namespace gqzip {

std::vector<uint8_t> CompressionEngine::pack_dna(const std::vector<std::string>& sequences) {
    std::vector<uint8_t> out;
    // Total sequences count (4 bytes)
    uint32_t count = static_cast<uint32_t>(sequences.size());
    out.resize(4);
    std::memcpy(out.data(), &count, 4);

    for (const auto& seq : sequences) {
        uint32_t len = static_cast<uint32_t>(seq.size());
        size_t offset = out.size();
        out.resize(offset + 4);
        std::memcpy(out.data() + offset, &len, 4);

        // 2-bit packing: 4 bases per byte
        size_t packed_bytes = (len + 3) / 4;
        std::vector<uint8_t> packed(packed_bytes, 0);

        std::vector<uint32_t> n_positions;
        for (size_t i = 0; i < len; ++i) {
            uint8_t code = 0;
            switch (seq[i]) {
                case 'A': case 'a': code = 0; break;
                case 'C': case 'c': code = 1; break;
                case 'G': case 'g': code = 2; break;
                case 'T': case 't': code = 3; break;
                default:            code = 0; n_positions.push_back(static_cast<uint32_t>(i)); break;
            }
            packed[i / 4] |= (code << ((i % 4) * 2));
        }

        out.insert(out.end(), packed.begin(), packed.end());

        // Store N positions (Patent Claim 2: Sparse N-mask vectorization)
        if (n_positions.empty()) {
            out.push_back(0x00);
        } else {
            out.push_back(0xFF);
            uint32_t num_n = static_cast<uint32_t>(n_positions.size());
            size_t n_off = out.size();
            out.resize(n_off + 4 + num_n * 4);
            std::memcpy(out.data() + n_off, &num_n, 4);
            std::memcpy(out.data() + n_off + 4, n_positions.data(), num_n * 4);
        }
    }
    return out;
}

std::vector<std::string> CompressionEngine::unpack_dna(const uint8_t* data, size_t size, size_t record_count) {
    std::vector<std::string> sequences;
    if (size < 4) return sequences;

    uint32_t count = 0;
    std::memcpy(&count, data, 4);
    sequences.reserve(count);

    size_t offset = 4;
    for (size_t r = 0; r < count && offset + 4 <= size; ++r) {
        uint32_t len = 0;
        std::memcpy(&len, data + offset, 4);
        offset += 4;

        size_t packed_bytes = (len + 3) / 4;
        if (offset + packed_bytes > size) break;

        const uint8_t* packed_ptr = data + offset;
        offset += packed_bytes;

        static const char BASES[4] = {'A', 'C', 'G', 'T'};
        std::string seq;
        seq.resize(len);

        for (size_t i = 0; i < len; ++i) {
            uint8_t code = (packed_ptr[i / 4] >> ((i % 4) * 2)) & 0x03;
            seq[i] = BASES[code];
        }

        if (offset < size) {
            uint8_t n_flag = data[offset++];
            if (n_flag == 0xFF && offset + 4 <= size) {
                uint32_t num_n = 0;
                std::memcpy(&num_n, data + offset, 4);
                offset += 4;

                for (size_t k = 0; k < num_n && offset + 4 <= size; ++k) {
                    uint32_t n_pos = 0;
                    std::memcpy(&n_pos, data + offset, 4);
                    offset += 4;
                    if (n_pos < len) {
                        seq[n_pos] = 'N';
                    }
                }
            }
        }

        sequences.push_back(std::move(seq));
    }
    return sequences;
}

std::vector<uint8_t> CompressionEngine::encode_headers(const std::vector<std::string>& headers) {
    if (headers.empty()) return {};

    // Patent Claim 2: Header Prefix-Delta Tokenization
    const std::string& first = headers[0];
    size_t common_len = first.size();
    for (size_t i = 1; i < headers.size(); ++i) {
        size_t j = 0;
        while (j < common_len && j < headers[i].size() && first[j] == headers[i][j]) {
            j++;
        }
        common_len = j;
    }

    std::string flat;
    if (common_len > 8) {
        flat += "P:";
        flat += first.substr(0, common_len);
        flat += '\n';
        for (const auto& h : headers) {
            flat += h.substr(common_len);
            flat += '\n';
        }
    } else {
        flat += "R:\n";
        for (const auto& h : headers) {
            flat += h;
            flat += '\n';
        }
    }

    return RansCodec::compress_stream(reinterpret_cast<const uint8_t*>(flat.data()), flat.size());
}

std::vector<std::string> CompressionEngine::decode_headers(const uint8_t* data, size_t size, size_t record_count) {
    auto decomp = RansCodec::decompress_stream(data, size, record_count * 64);
    std::string_view flat(reinterpret_cast<const char*>(decomp.data()), decomp.size());

    std::vector<std::string> headers;
    headers.reserve(record_count);

    size_t start = 0;
    std::string_view prefix;
    
    size_t first_nl = flat.find('\n');
    if (first_nl != std::string_view::npos && flat.starts_with("P:")) {
        prefix = flat.substr(2, first_nl - 2);
        start = first_nl + 1;
    } else if (first_nl != std::string_view::npos && flat.starts_with("R:\n")) {
        start = first_nl + 1;
    }

    while (start < flat.size() && headers.size() < record_count) {
        size_t next_nl = flat.find('\n', start);
        if (next_nl == std::string_view::npos) {
            std::string h(prefix);
            h += flat.substr(start);
            headers.push_back(std::move(h));
            break;
        }
        std::string h(prefix);
        h += flat.substr(start, next_nl - start);
        headers.push_back(std::move(h));
        start = next_nl + 1;
    }
    return headers;
}

BlockData CompressionEngine::compress_batch(const FastqBatch& batch) {
    BlockData block;
    block.record_count = static_cast<uint32_t>(batch.size());
    block.binning_level = m_options.binning_level;
    block.flags = (m_options.lossless || m_options.binning_level == BinningLevel::LEVEL_5_LOSSLESS) ? 0x0001 : 0x0000;

    // Headers
    block.header_stream = encode_headers(batch.headers);

    // DNA
    block.dna_stream = pack_dna(batch.sequences);

    // Quality Scores & Residuals
    auto quant_results = m_quantizer.quantize_batch(batch);
    std::string flat_quality;
    std::vector<int8_t> flat_residuals;

    for (const auto& qr : quant_results) {
        flat_quality += qr.quantized_quality;
        flat_quality += '\n';
        if (block.flags & 0x0001) {
            flat_residuals.insert(flat_residuals.end(), qr.residuals.begin(), qr.residuals.end());
        }
    }

    block.quality_stream = RansCodec::compress_stream(
        reinterpret_cast<const uint8_t*>(flat_quality.data()),
        flat_quality.size()
    );

    if (block.flags & 0x0001) {
        block.residual_stream = RansCodec::compress_residuals(flat_residuals);
    }

    // Compute checksum
    uint32_t crc = CRC32::calculate(block.header_stream.data(), block.header_stream.size());
    crc = CRC32::calculate(block.dna_stream.data(), block.dna_stream.size(), crc);
    crc = CRC32::calculate(block.quality_stream.data(), block.quality_stream.size(), crc);
    if (!block.residual_stream.empty()) {
        crc = CRC32::calculate(block.residual_stream.data(), block.residual_stream.size(), crc);
    }
    block.crc32 = crc;

    return block;
}

FastqBatch CompressionEngine::decompress_block(const BlockData& block) {
    FastqBatch batch;
    batch.reserve(block.record_count);

    auto headers = decode_headers(block.header_stream.data(), block.header_stream.size(), block.record_count);
    auto sequences = unpack_dna(block.dna_stream.data(), block.dna_stream.size(), block.record_count);

    auto decomp_qual = RansCodec::decompress_stream(
        block.quality_stream.data(),
        block.quality_stream.size(),
        block.record_count * 150
    );
    std::string_view flat_qual(reinterpret_cast<const char*>(decomp_qual.data()), decomp_qual.size());

    std::vector<std::string_view> qual_lines;
    size_t start = 0;
    while (start < flat_qual.size() && qual_lines.size() < block.record_count) {
        size_t next_nl = flat_qual.find('\n', start);
        if (next_nl == std::string_view::npos) {
            qual_lines.push_back(flat_qual.substr(start));
            break;
        }
        qual_lines.push_back(flat_qual.substr(start, next_nl - start));
        start = next_nl + 1;
    }

    std::vector<int8_t> all_residuals;
    if (block.flags & 0x0001 && !block.residual_stream.empty()) {
        all_residuals = RansCodec::decompress_residuals(
            block.residual_stream.data(),
            block.residual_stream.size(),
            block.record_count * 150
        );
    }

    size_t res_offset = 0;
    for (size_t i = 0; i < block.record_count; ++i) {
        std::string h = (i < headers.size()) ? headers[i] : ("@READ_" + std::to_string(i));
        std::string s = (i < sequences.size()) ? sequences[i] : "";
        std::string q;

        if (i < qual_lines.size()) {
            std::string_view q_line = qual_lines[i];
            if (block.flags & 0x0001 && res_offset + q_line.size() <= all_residuals.size()) {
                std::vector<int8_t> sub_res(
                    all_residuals.begin() + res_offset,
                    all_residuals.begin() + res_offset + q_line.size()
                );
                q = Quantizer::restore(q_line, sub_res);
                res_offset += q_line.size();
            } else {
                q = std::string(q_line);
            }
        }

        batch.add(std::move(h), std::move(s), std::move(q));
    }

    return batch;
}

CompressionStats CompressionEngine::compress(std::istream& in, std::ostream& out) {
    auto t_start = std::chrono::high_resolution_clock::now();
    CompressionStats stats;

    // Write Global Magic & Header
    out.write(reinterpret_cast<const char*>(&MAGIC_HEADER), 4);
    out.write(reinterpret_cast<const char*>(&FORMAT_VERSION), 2);
    uint8_t b_level = static_cast<uint8_t>(m_options.binning_level);
    out.write(reinterpret_cast<const char*>(&b_level), 1);
    uint8_t flags = m_options.lossless ? 1 : 0;
    out.write(reinterpret_cast<const char*>(&flags), 1);

    FastqBatch batch;
    while (SIMDScanner::read_batch(in, batch, m_options.block_records)) {
        for (size_t i = 0; i < batch.size(); ++i) {
            stats.raw_fastq_bytes += batch.headers[i].size() + batch.sequences[i].size() + batch.qualities[i].size() + 4;
            stats.raw_quality_bytes += batch.qualities[i].size();
        }
        stats.total_records += batch.size();

        BlockData block = compress_batch(batch);
        block.serialize(out);

        stats.compressed_quality_bytes += block.quality_stream.size() + block.residual_stream.size();
        stats.compressed_bytes += block.header_stream.size() + block.dna_stream.size() + block.quality_stream.size() + block.residual_stream.size() + 32;
    }

    // Write End-of-Stream Marker (Block with 0 records)
    BlockData eos;
    eos.record_count = 0;
    eos.serialize(out);

    auto t_end = std::chrono::high_resolution_clock::now();
    stats.compression_time_sec = std::chrono::duration<double>(t_end - t_start).count();
    return stats;
}

CompressionStats CompressionEngine::decompress(std::istream& in, std::ostream& out) {
    auto t_start = std::chrono::high_resolution_clock::now();
    CompressionStats stats;

    uint32_t magic = 0;
    if (!in.read(reinterpret_cast<char*>(&magic), 4) || magic != MAGIC_HEADER) {
        throw std::runtime_error("Invalid file format: missing GQZ magic header");
    }

    uint16_t version = 0;
    in.read(reinterpret_cast<char*>(&version), 2);
    uint8_t b_level = 0, flags = 0;
    in.read(reinterpret_cast<char*>(&b_level), 1);
    in.read(reinterpret_cast<char*>(&flags), 1);

    BlockData block;
    while (block.deserialize(in) && block.record_count > 0) {
        FastqBatch batch = decompress_block(block);
        for (size_t i = 0; i < batch.size(); ++i) {
            out << batch.headers[i] << "\n"
                << batch.sequences[i] << "\n"
                << "+\n"
                << batch.qualities[i] << "\n";
            stats.raw_fastq_bytes += batch.headers[i].size() + batch.sequences[i].size() + batch.qualities[i].size() + 4;
        }
        stats.total_records += batch.size();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    stats.decompression_time_sec = std::chrono::duration<double>(t_end - t_start).count();
    return stats;
}

} // namespace gqzip
