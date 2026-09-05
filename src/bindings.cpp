#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gqzip/engine.hpp"
#include "gqzip/rans_codec.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_gqzip_cpp, m) {
    m.doc() = "GQZip High-Performance C++ Compression Extension";

    py::enum_<gqzip::BinningLevel>(m, "BinningLevel")
        .value("LEVEL_1_ILLUMINA_8BIN", gqzip::BinningLevel::LEVEL_1_ILLUMINA_8BIN)
        .value("LEVEL_2_COARSE_4BIN", gqzip::BinningLevel::LEVEL_2_COARSE_4BIN)
        .value("LEVEL_3_ADAPTIVE_CONTEXT", gqzip::BinningLevel::LEVEL_3_ADAPTIVE_CONTEXT)
        .value("LEVEL_4_BINARY", gqzip::BinningLevel::LEVEL_4_BINARY)
        .value("LEVEL_5_LOSSLESS", gqzip::BinningLevel::LEVEL_5_LOSSLESS)
        .export_values();

    py::class_<gqzip::CompressionOptions>(m, "CompressionOptions")
        .def(py::init<>())
        .def_readwrite("binning_level", &gqzip::CompressionOptions::binning_level)
        .def_readwrite("lossless", &gqzip::CompressionOptions::lossless)
        .def_readwrite("num_threads", &gqzip::CompressionOptions::num_threads)
        .def_readwrite("verbose", &gqzip::CompressionOptions::verbose);

    py::class_<gqzip::CompressionStats>(m, "CompressionStats")
        .def(py::init<>())
        .def_readwrite("total_records", &gqzip::CompressionStats::total_records)
        .def_readwrite("raw_fastq_bytes", &gqzip::CompressionStats::raw_fastq_bytes)
        .def_readwrite("compressed_bytes", &gqzip::CompressionStats::compressed_bytes)
        .def_readwrite("compression_time_sec", &gqzip::CompressionStats::compression_time_sec)
        .def("compression_ratio", &gqzip::CompressionStats::compression_ratio)
        .def("throughput_mb_s", &gqzip::CompressionStats::throughput_mb_s);

    py::class_<gqzip::RansCodec>(m, "RansCodec")
        .def(py::init<>())
        .def_static("encode", [](const std::string& data) {
            auto vec = gqzip::RansCodec::encode(reinterpret_cast<const uint8_t*>(data.data()), data.size());
            return py::bytes(reinterpret_cast<const char*>(vec.data()), vec.size());
        })
        .def_static("decode", [](const std::string& data, size_t orig_sz) {
            auto vec = gqzip::RansCodec::decode(reinterpret_cast<const uint8_t*>(data.data()), data.size(), orig_sz);
            return py::bytes(reinterpret_cast<const char*>(vec.data()), vec.size());
        })
        .def_static("encode_order1", [](const std::string& data) {
            auto vec = gqzip::RansCodec::encode_order1(reinterpret_cast<const uint8_t*>(data.data()), data.size());
            return py::bytes(reinterpret_cast<const char*>(vec.data()), vec.size());
        })
        .def_static("decode_order1", [](const std::string& data, size_t orig_sz) {
            auto vec = gqzip::RansCodec::decode_order1(reinterpret_cast<const uint8_t*>(data.data()), data.size(), orig_sz);
            return py::bytes(reinterpret_cast<const char*>(vec.data()), vec.size());
        });
}
