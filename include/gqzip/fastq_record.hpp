#pragma once

#include "common.hpp"
#include <string>
#include <vector>
#include <span>

namespace gqzip {

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Columnar batch representation for maximum SIMD efficiency and cache locality
struct FastqBatch {
    std::vector<std::string> headers;
    std::vector<std::string> sequences;
    std::vector<std::string> qualities;

    size_t size() const {
        return headers.size();
    }

    bool empty() const {
        return headers.empty();
    }

    void reserve(size_t n) {
        headers.reserve(n);
        sequences.reserve(n);
        qualities.reserve(n);
    }

    void clear() {
        headers.clear();
        sequences.clear();
        qualities.clear();
    }

    void add(std::string h, std::string s, std::string q) {
        headers.push_back(std::move(h));
        sequences.push_back(std::move(s));
        qualities.push_back(std::move(q));
    }
};

} // namespace gqzip
