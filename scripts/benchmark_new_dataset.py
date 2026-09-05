#!/usr/bin/env python3
import os
import sys
import time
import hashlib

# Ensure python directory is in sys.path
sys.path.insert(0, os.path.abspath("python"))
from gqzip.engine import CompressionEngine
from gqzip.core import CompressionOptions, BinningLevel

input_path = r"C:\Users\jtaylor1\Downloads\ERR3239334_clean_records.fastq"
raw_bytes = os.path.getsize(input_path)

with open(input_path, "rb") as f:
    raw_hash = hashlib.sha256(f.read()).hexdigest().upper()

print("==========================================================================")
print("     BENCHMARKING GQZIP ACROSS ALL 3 MODES ON NEW DATASET (ERR3239334)")
print("==========================================================================")
print(f"Input File: {input_path}")
print(f"Raw Input Size: {raw_bytes:,} bytes ({raw_bytes / (1024*1024):.2f} MB)")
print(f"Raw Input SHA-256: {raw_hash}\n")

modes = [
    (5, "Mode -b 5 (Lossless Reversible)", True),
    (3, "Mode -b 3 (Adaptive Context)", False),
    (4, "Mode -b 4 (Binary 2-Bin)", False),
]

results = []

for b_level, name, is_lossless in modes:
    out_path = f"scratch/ERR3239334_mode{b_level}.gqz"
    os.makedirs("scratch", exist_ok=True)
    
    options = CompressionOptions(
        binning_level=BinningLevel(b_level),
        lossless=is_lossless,
        verbose=False
    )
    engine = CompressionEngine(options)
    
    t0 = time.time()
    with open(input_path, "rb") as in_f, open(out_path, "wb") as out_f:
        stats = engine.compress_stream(in_f, out_f)
    t1 = time.time()
    
    comp_bytes = os.path.getsize(out_path)
    ratio = raw_bytes / comp_bytes if comp_bytes > 0 else 1.0
    speed_mb_s = (raw_bytes / (1024 * 1024)) / (t1 - t0)
    
    with open(out_path, "rb") as f:
        gqz_hash = hashlib.sha256(f.read()).hexdigest().upper()
        
    # If Lossless (-b 5), test roundtrip decompression
    decomp_hash = "N/A (Lossy Mode)"
    is_bit_exact = "N/A"
    if is_lossless:
        restored_path = "scratch/ERR3239334_restored_mode5.fastq"
        with open(out_path, "rb") as in_f, open(restored_path, "wb") as out_f:
            engine.decompress_stream(in_f, out_f)
        with open(restored_path, "rb") as f:
            decomp_hash = hashlib.sha256(f.read()).hexdigest().upper()
        is_bit_exact = "100.00% PASS (BIT-EXACT MATCH)" if decomp_hash == raw_hash else "FAIL"

    results.append({
        "name": name,
        "mode": b_level,
        "comp_bytes": comp_bytes,
        "ratio": ratio,
        "speed": speed_mb_s,
        "time": t1 - t0,
        "gqz_hash": gqz_hash,
        "decomp_hash": decomp_hash,
        "is_bit_exact": is_bit_exact
    })

print("==========================================================================")
print("                        COMPRESSION BENCHMARK RESULTS")
print("==========================================================================")
for r in results:
    print(f"\n--- {r['name']} ---")
    print(f"  Compressed Size : {r['comp_bytes']:,} bytes ({r['comp_bytes']/(1024*1024):.2f} MB)")
    print(f"  Compression Ratio: {r['ratio']:.2f}x ({100*(1-r['comp_bytes']/raw_bytes):.2f}% space savings)")
    print(f"  Encode Speed     : {r['speed']:.2f} MB/s (Time: {r['time']:.3f}s)")
    print(f"  Container SHA-256: {r['gqz_hash']}")
    if r["mode"] == 5:
        print(f"  Restored SHA-256 : {r['decomp_hash']}")
        print(f"  Bit-Exact Status : {r['is_bit_exact']}")
print("==========================================================================")
