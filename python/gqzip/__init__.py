"""
gqzip: High-Performance, Context-Aware FASTQ Quality Score Compression Engine.
"""

from .core import BinningLevel, CompressionOptions, CompressionStats, FASTQRecord
from .quantizer import Quantizer, calculate_local_entropy, calculate_cycle_decay, compute_lloyd_max_centroids
from .rans_codec import RansCodec
from .engine import CompressionEngine, compress_file, decompress_file
from .s3_stream import S3StreamingGateway, compress_to_s3, decompress_from_s3

__version__ = "1.1.0"
__all__ = [
    "BinningLevel",
    "CompressionOptions",
    "CompressionStats",
    "FASTQRecord",
    "Quantizer",
    "RansCodec",
    "CompressionEngine",
    "S3StreamingGateway",
    "compress_file",
    "decompress_file",
    "compress_to_s3",
    "decompress_from_s3",
    "calculate_local_entropy",
    "calculate_cycle_decay",
    "compute_lloyd_max_centroids",
]
