#include "gqzip/common.hpp"
#include "gqzip/engine.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#endif

void print_help(const char* prog) {
    std::cout << "gqzip: High-Performance Context-Aware Quality Score Compression Engine\n"
              << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  -c, --compress         Compress FASTQ to .gqz (default)\n"
              << "  -d, --decompress       Decompress .gqz to standard FASTQ\n"
              << "  -b, --binning <1-5>    Binning level (1=Illumina 8-bin, 2=Coarse 4-bin,\n"
              << "                         3=Context-Adaptive (Default), 4=Binary 2-bin,\n"
              << "                         5=Lossless)\n"
              << "  -l, --lossless         Enable 100% bit-exact lossless residual restoration\n"
              << "  -t, --threads <N>      Number of worker threads (default: 4)\n"
              << "  -i, --input <file>     Input file (default: stdin)\n"
              << "  -o, --output <file>    Output file (default: stdout)\n"
              << "  -v, --verbose          Print compression statistics and throughput\n"
              << "  -h, --help             Display this help message and exit\n\n"
              << "Examples:\n"
              << "  " << prog << " -c -b 3 -i input.fastq -o output.gqz\n"
              << "  " << prog << " -d -i output.gqz -o reconstructed.fastq\n"
              << "  cat reads.fastq | " << prog << " -c -b 3 | " << prog << " -d > restored.fastq\n";
}

int main(int argc, char* argv[]) {
    // Set binary mode for stdin/stdout on Windows
#if defined(_WIN32)
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    gqzip::CompressionOptions options;
    bool decompress_mode = false;
    std::string input_path;
    std::string output_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help(argv[0]);
            return 0;
        } else if (arg == "-c" || arg == "--compress") {
            decompress_mode = false;
        } else if (arg == "-d" || arg == "--decompress") {
            decompress_mode = true;
        } else if (arg == "-l" || arg == "--lossless") {
            options.lossless = true;
            options.binning_level = gqzip::BinningLevel::LEVEL_5_LOSSLESS;
        } else if (arg == "-v" || arg == "--verbose") {
            options.verbose = true;
        } else if ((arg == "-b" || arg == "--binning") && i + 1 < argc) {
            int lvl = std::stoi(argv[++i]);
            if (lvl >= 1 && lvl <= 5) {
                options.binning_level = static_cast<gqzip::BinningLevel>(lvl);
                if (lvl == 5) options.lossless = true;
            }
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            options.num_threads = std::stoul(argv[++i]);
        } else if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            input_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        }
    }

    // Input stream setup
    std::ifstream file_in;
    std::istream* in_ptr = &std::cin;
    if (!input_path.empty()) {
        file_in.open(input_path, std::ios::binary);
        if (!file_in.is_open()) {
            std::cerr << "Error: Could not open input file: " << input_path << "\n";
            return 1;
        }
        in_ptr = &file_in;
    }

    // Output stream setup
    std::ofstream file_out;
    std::ostream* out_ptr = &std::cout;
    if (!output_path.empty()) {
        file_out.open(output_path, std::ios::binary);
        if (!file_out.is_open()) {
            std::cerr << "Error: Could not open output file: " << output_path << "\n";
            return 1;
        }
        out_ptr = &file_out;
    }

    gqzip::CompressionEngine engine(options);

    try {
        if (!decompress_mode) {
            auto stats = engine.compress(*in_ptr, *out_ptr);
            if (options.verbose) {
                std::cerr << "=== gqzip Compression Complete ===\n"
                          << "Total Reads Processed   : " << stats.total_records << "\n"
                          << "Raw FASTQ Size          : " << stats.raw_fastq_bytes / (1024.0 * 1024.0) << " MB\n"
                          << "Compressed Size         : " << stats.compressed_bytes / (1024.0 * 1024.0) << " MB\n"
                          << "Overall Ratio           : " << stats.compression_ratio() << "x\n"
                          << "Quality Footprint Saving: " << stats.quality_space_saving_pct() << "%\n"
                          << "Throughput              : " << stats.throughput_mb_s() << " MB/s\n";
            }
        } else {
            auto stats = engine.decompress(*in_ptr, *out_ptr);
            if (options.verbose) {
                std::cerr << "=== gqzip Decompression Complete ===\n"
                          << "Total Reads Restored: " << stats.total_records << "\n"
                          << "Decompression Time  : " << stats.decompression_time_sec << " s\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "gqzip Engine Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
