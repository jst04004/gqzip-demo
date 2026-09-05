#!/usr/bin/env python3
"""
Comprehensive Local Stress, Edge-Case, and Fuzz Test Suite for GQZip.
Tests error handling, diverse headers, long reads, corrupted input, and memory bounds.
"""

import os
import sys
import tempfile
import hashlib
import time
import random

sys.path.insert(0, os.path.abspath("python"))
from gqzip.engine import CompressionEngine
from gqzip.core import CompressionOptions, BinningLevel, FASTQRecord

print("==========================================================================")
print("       GQZIP COMPREHENSIVE LOCAL STRESS, FUZZ & EDGE-CASE SUITE")
print("==========================================================================")

test_dir = tempfile.mkdtemp(prefix="gqzip_stress_suite_")
print(f"Working Test Directory: {test_dir}\n")

test_results = []

def record_test(name, status, detail=""):
    test_results.append({"name": name, "status": status, "detail": detail})
    symbol = "[PASS]" if status == "PASS" else "[FAIL]"
    print(f"{symbol} {name}")
    if detail:
        print(f"         Detail: {detail}")

# ------------------------------------------------------------------------------
# Test 1: Empty & Zero-Byte FASTQ Files
# ------------------------------------------------------------------------------
try:
    empty_in = os.path.join(test_dir, "empty.fastq")
    empty_out = os.path.join(test_dir, "empty.gqz")
    empty_rest = os.path.join(test_dir, "empty_restored.fastq")
    
    with open(empty_in, "wb") as f:
        pass
        
    engine = CompressionEngine(CompressionOptions(lossless=True))
    with open(empty_in, "rb") as in_f, open(empty_out, "wb") as out_f:
        stats = engine.compress_stream(in_f, out_f)
        
    with open(empty_out, "rb") as in_f, open(empty_rest, "wb") as out_f:
        engine.decompress_stream(in_f, out_f)
        
    if os.path.getsize(empty_rest) == 0:
        record_test("Test 1: Empty (0-Byte) File Handling", "PASS", "Handled 0-byte file without crashing.")
    else:
        record_test("Test 1: Empty (0-Byte) File Handling", "FAIL", f"Expected 0-byte output, got {os.path.getsize(empty_rest)} bytes.")
except Exception as e:
    record_test("Test 1: Empty (0-Byte) File Handling", "FAIL", str(e))

# ------------------------------------------------------------------------------
# Test 2: Ultra-Long Read Lengths (10,000 bp Long-Read PacBio Simulation)
# ------------------------------------------------------------------------------
try:
    long_in = os.path.join(test_dir, "long_reads.fastq")
    long_out = os.path.join(test_dir, "long_reads.gqz")
    long_rest = os.path.join(test_dir, "long_reads_restored.fastq")
    
    bases = ["A", "C", "G", "T"]
    long_records = []
    with open(long_in, "w", newline="\n") as f:
        for r_id in range(10):
            seq = "".join(random.choices(bases, k=10000))
            qual = "".join(chr(random.randint(33, 73)) for _ in range(10000))
            f.write(f"@LONG_READ_{r_id+1}_LEN_10000\n{seq}\n+\n{qual}\n")
            
    with open(long_in, "rb") as f:
        orig_hash = hashlib.sha256(f.read()).hexdigest().upper()
        
    engine = CompressionEngine(CompressionOptions(lossless=True, binning_level=BinningLevel.LEVEL_5_LOSSLESS))
    with open(long_in, "rb") as in_f, open(long_out, "wb") as out_f:
        engine.compress_stream(in_f, out_f)
        
    with open(long_out, "rb") as in_f, open(long_rest, "wb") as out_f:
        engine.decompress_stream(in_f, out_f)
        
    with open(long_rest, "rb") as f:
        rest_hash = hashlib.sha256(f.read()).hexdigest().upper()
        
    if orig_hash == rest_hash:
        record_test("Test 2: Ultra-Long Reads (10,000 bp)", "PASS", "100.00% SHA-256 bit-exact match on 10,000 bp reads.")
    else:
        record_test("Test 2: Ultra-Long Reads (10,000 bp)", "FAIL", "Checksum mismatch on 10,000 bp reads.")
except Exception as e:
    record_test("Test 2: Ultra-Long Reads (10,000 bp)", "FAIL", str(e))

# ------------------------------------------------------------------------------
# Test 3: Non-Standard Arbitrary Header Formats (Element / MGI / Custom Text)
# ------------------------------------------------------------------------------
try:
    cust_in = os.path.join(test_dir, "custom_headers.fastq")
    cust_out = os.path.join(test_dir, "custom_headers.gqz")
    cust_rest = os.path.join(test_dir, "custom_headers_restored.fastq")
    
    headers = [
        "@ELEMENT_AVITI_RUN_99_FLOWCELL_A1:LANE1:READ1_X192_Y304",
        "@MGI_SEQ_DNB_2026_COL_001_ROW_4022",
        "@CUSTOM_RESEARCH_SAMPLE_ID_99824_BARCODE_ATCGATCG",
        "@SIMPLE_ID_1"
    ]
    with open(cust_in, "w", newline="\n") as f:
        for idx, h in enumerate(headers):
            f.write(f"{h}\nACGTACGTACGT\n+\n????????????\n")
            
    with open(cust_in, "rb") as f:
        orig_hash = hashlib.sha256(f.read()).hexdigest().upper()
        
    engine = CompressionEngine(CompressionOptions(lossless=True))
    with open(cust_in, "rb") as in_f, open(cust_out, "wb") as out_f:
        engine.compress_stream(in_f, out_f)
        
    with open(cust_out, "rb") as in_f, open(cust_rest, "wb") as out_f:
        engine.decompress_stream(in_f, out_f)
        
    with open(cust_rest, "rb") as f:
        rest_hash = hashlib.sha256(f.read()).hexdigest().upper()
        
    if orig_hash == rest_hash:
        record_test("Test 3: Non-Standard Headers (MGI / Element)", "PASS", "100.00% SHA-256 match across non-Illumina header styles.")
    else:
        record_test("Test 3: Non-Standard Headers (MGI / Element)", "FAIL", "Checksum mismatch on non-standard headers.")
except Exception as e:
    record_test("Test 3: Non-Standard Headers (MGI / Element)", "FAIL", str(e))

# ------------------------------------------------------------------------------
# Test 4: Corrupted Container / Fuzz Injection Protection
# ------------------------------------------------------------------------------
try:
    corrupt_out = os.path.join(test_dir, "corrupted.gqz")
    corrupt_rest = os.path.join(test_dir, "corrupted_restored.fastq")
    
    # Write invalid corrupted binary payload
    with open(corrupt_out, "wb") as f:
        f.write(b"GQZ\x01INVALID_CORRUPTED_BYTES_HERE_1234567890")
        
    engine = CompressionEngine()
    caught = False
    try:
        with open(corrupt_out, "rb") as in_f, open(corrupt_rest, "wb") as out_f:
            engine.decompress_stream(in_f, out_f)
    except (ValueError, Exception):
        caught = True
        
    if caught:
        record_test("Test 4: Fuzz Injection & Corrupted Container Check", "PASS", "Caught corrupted payload safely without infinite loops or silent failure.")
    else:
        record_test("Test 4: Fuzz Injection & Corrupted Container Check", "FAIL", "Engine did not detect corrupted payload.")
except Exception as e:
    record_test("Test 4: Fuzz Injection & Corrupted Container Check", "FAIL", str(e))

# ------------------------------------------------------------------------------
# Test 5: Memory Boundedness Check (<50 MB RSS Ceiling)
# ------------------------------------------------------------------------------
try:
    mem_in = os.path.join(test_dir, "mem_test.fastq")
    mem_out = os.path.join(test_dir, "mem_test.gqz")
    
    # Generate 20,000 FASTQ records
    with open(mem_in, "w", newline="\n") as f:
        for idx in range(20000):
            f.write(f"@READ_{idx}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\n????????????????????????????????\n")
            
    engine = CompressionEngine(CompressionOptions(binning_level=BinningLevel.LEVEL_3_ADAPTIVE_CONTEXT))
    with open(mem_in, "rb") as in_f, open(mem_out, "wb") as out_f:
        stats = engine.compress_stream(in_f, out_f)
        
    if stats.peak_memory_mb < 50.0:
        record_test("Test 5: Memory RSS Bound Check (<50 MB)", "PASS", f"Peak RSS was {stats.peak_memory_mb:.1f} MB (well below 50 MB ceiling).")
    else:
        record_test("Test 5: Memory RSS Bound Check (<50 MB)", "FAIL", f"Peak RSS was {stats.peak_memory_mb:.1f} MB.")
except Exception as e:
    record_test("Test 5: Memory RSS Bound Check (<50 MB)", "FAIL", str(e))

print("\n==========================================================================")
pass_count = sum(1 for r in test_results if r["status"] == "PASS")
total_count = len(test_results)
print(f"SUITE VERDICT: Passed {pass_count}/{total_count} Tests ({100.0 * pass_count / total_count:.1f}%)")
print("==========================================================================")
