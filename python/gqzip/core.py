from enum import IntEnum
from dataclasses import dataclass
from typing import List, Optional

MAGIC_HEADER = b"GQZ\x01"
FORMAT_VERSION = 1
DEFAULT_BLOCK_RECORDS = 50000

class BinningLevel(IntEnum):
    LEVEL_1_ILLUMINA8 = 1       # Standard Illumina 8-bin scheme
    LEVEL_2_COARSE4 = 2         # Coarse 4-bin scheme
    LEVEL_3_ADAPTIVE_CONTEXT = 3 # Dynamic context-aware 1D window quantizer (Default)
    LEVEL_4_BINARY = 4          # Binary high/low quality thresholding
    LEVEL_5_LOSSLESS = 5        # 100% Bit-exact lossless restoration with residual stream

@dataclass
class FASTQRecord:
    header: str
    sequence: str
    quality: str

@dataclass
class CompressionOptions:
    binning_level: BinningLevel = BinningLevel.LEVEL_3_ADAPTIVE_CONTEXT
    lossless: bool = False
    block_records: int = DEFAULT_BLOCK_RECORDS
    num_threads: int = 4
    verbose: bool = False

@dataclass
class CompressionStats:
    total_records: int = 0
    raw_fastq_bytes: int = 0
    compressed_bytes: int = 0
    raw_quality_bytes: int = 0
    compressed_quality_bytes: int = 0
    compression_time_sec: float = 0.0
    decompression_time_sec: float = 0.0
    peak_memory_mb: float = 0.0

    @property
    def compression_ratio(self) -> float:
        return (self.raw_fastq_bytes / self.compressed_bytes) if self.compressed_bytes > 0 else 1.0

    @property
    def space_saving_pct(self) -> float:
        if self.raw_fastq_bytes <= 0:
            return 0.0
        return (1.0 - self.compressed_bytes / self.raw_fastq_bytes) * 100.0

    @property
    def quality_saving_pct(self) -> float:
        if self.raw_quality_bytes <= 0:
            return 0.0
        return (1.0 - self.compressed_quality_bytes / self.raw_quality_bytes) * 100.0

    @property
    def throughput_mb_s(self) -> float:
        if self.compression_time_sec <= 0:
            return 0.0
        return (self.raw_fastq_bytes / (1024.0 * 1024.0)) / self.compression_time_sec
