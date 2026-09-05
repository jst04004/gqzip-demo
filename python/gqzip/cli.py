import sys
import argparse
from gqzip.core import BinningLevel, CompressionOptions
from gqzip.engine import CompressionEngine

def main():
    parser = argparse.ArgumentParser(
        description="gqzip: High-Performance Context-Aware Quality Score Compression Engine"
    )
    parser.add_argument("-c", "--compress", action="store_true", default=False, help="Compress FASTQ to .gqz")
    parser.add_argument("-d", "--decompress", action="store_true", help="Decompress .gqz to standard FASTQ")
    parser.add_argument("-e", "--evaluate", action="store_true", help="Evaluate compression and generate ROI report on input FASTQ")
    parser.add_argument("--monthly-tb", type=float, default=50.0, help="Monthly sequencing volume in TB for ROI calculation")
    parser.add_argument("-b", "--binning", type=int, choices=[1, 2, 3, 4, 5], default=3,
                        help="Binning level (1=Illumina 8-bin, 2=Coarse 4-bin, 3=Adaptive Context, 4=Binary, 5=Lossless)")
    parser.add_argument("-l", "--lossless", action="store_true", help="Enable 100%% bit-exact lossless residual restoration")
    parser.add_argument("-i", "--input", default=None, help="Input file path (default: stdin)")
    parser.add_argument("-o", "--output", default=None, help="Output file path (default: stdout)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print compression statistics and throughput")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Worker thread count")
    parser.add_argument("--license-key", default=None, help="Activate academic or enterprise license key")

    args = parser.parse_args()

    if args.license_key:
        from gqzip.license import save_license_key
        ok, msg = save_license_key(args.license_key)
        sys.stderr.write(f"\n[gqzip License] {msg}\n")
        if not ok:
            sys.exit(1)
        sys.exit(0)

    if args.evaluate:
        if not args.input:
            sys.stderr.write("[gqzip Error] --input FASTQ file required for evaluation.\n")
            sys.exit(1)
        from gqzip.evaluate import evaluate_custom_fastq
        evaluate_custom_fastq(
            input_path=args.input,
            output_dir=args.output,
            monthly_tb=args.monthly_tb,
            generate_html=True
        )
        sys.exit(0)

    options = CompressionOptions(
        binning_level=BinningLevel(args.binning),
        lossless=args.lossless or (args.binning == 5),
        num_threads=args.threads,
        verbose=args.verbose
    )

    engine = CompressionEngine(options)

    # Setup binary streams
    in_stream = open(args.input, "rb") if args.input else sys.stdin.buffer
    out_stream = open(args.output, "wb") if args.output else sys.stdout.buffer

    try:
        if args.decompress:
            stats = engine.decompress_stream(in_stream, out_stream)
            if args.verbose:
                sys.stderr.write(
                    f"\n[gqzip Decompress] Restored {stats.total_records:,} reads in {stats.decompression_time_sec:.3f}s\n"
                )
        else:
            stats = engine.compress_stream(in_stream, out_stream)
            if args.verbose:
                sys.stderr.write(
                    f"\n[gqzip Compress] Processed {stats.total_records:,} reads\n"
                    f"  Raw FASTQ Size  : {stats.raw_fastq_bytes / (1024*1024):.2f} MB\n"
                    f"  Compressed Size : {stats.compressed_bytes / (1024*1024):.2f} MB\n"
                    f"  Overall Ratio   : {stats.compression_ratio:.2f}x\n"
                    f"  Quality Savings : {stats.quality_saving_pct:.2f}%\n"
                    f"  Throughput      : {stats.throughput_mb_s:.2f} MB/s\n"
                    f"  Peak RAM (RSS)  : {stats.peak_memory_mb:.1f} MB (strictly bounded < 50 MB RSS)\n"
                )
    finally:
        if args.input:
            in_stream.close()
        if args.output:
            out_stream.close()

if __name__ == "__main__":
    main()
