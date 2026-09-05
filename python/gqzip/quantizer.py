import math
from typing import Tuple, List
from .core import BinningLevel

def calculate_local_entropy(sequence: str, pos: int, window_half_width: int = 3) -> float:
    """
    Computes local Shannon entropy H(W_i) over sliding window [pos - w, pos + w].
    H(W_i) = - sum(p * log2(p)) for p > 0.
    """
    seq_len = len(sequence)
    if seq_len == 0:
        return 0.0
        
    start = max(0, pos - window_half_width)
    end = min(seq_len, pos + window_half_width + 1)
    win_len = end - start
    if win_len == 0:
        return 0.0
        
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for k in range(start, end):
        base = sequence[k].upper()
        if base in counts:
            counts[base] += 1
            
    entropy = 0.0
    for c, count in counts.items():
        if count > 0:
            p = count / float(win_len)
            entropy -= p * math.log2(p)
    return entropy

def calculate_cycle_decay(
    pos: int,
    read_len: int,
    alpha: float = 0.35,
    gamma: float = 2.0,
    is_read_2: bool = False
) -> float:
    """Calculates empirical 3' cycle decay score D(i) with Paired-End Read 1/2 asymmetry."""
    if read_len <= 0:
        return 1.0
    # Read 2 in paired-end SBS exhibits ~35% faster optical decay and phasing accumulation
    effective_alpha = alpha * 1.35 if is_read_2 else alpha
    norm_pos = min(1.0, max(0.0, float(pos) / float(read_len)))
    decay = 1.0 - effective_alpha * math.pow(norm_pos, gamma)
    return max(0.0, min(1.0, decay))

def compute_lloyd_max_centroids(phred_counts: dict, num_bins: int = 4, max_iter: int = 20) -> List[int]:
    """Computes mathematically optimal 1D K-Means / Lloyd-Max centroids for empirical Phred distribution."""
    total_scores = sum(phred_counts.values())
    if total_scores == 0 or not phred_counts:
        return [6, 18, 30, 38] if num_bins == 4 else [10, 35]
        
    unique_phreds = sorted(phred_counts.keys())
    if len(unique_phreds) <= num_bins:
        centroids = sorted(unique_phreds)
        while len(centroids) < num_bins:
            centroids.append(centroids[-1] if centroids else 30)
        return centroids
        
    # Initialize centroids evenly across empirical range
    min_q = unique_phreds[0]
    max_q = unique_phreds[-1]
    step = (max_q - min_q) / float(num_bins + 1)
    centroids = [float(min_q + step * (i + 1)) for i in range(num_bins)]
    
    for _ in range(max_iter):
        # 1. Partition Phreds to nearest centroid
        clusters = [[] for _ in range(num_bins)]
        for q, count in phred_counts.items():
            nearest_idx = min(range(num_bins), key=lambda idx: abs(q - centroids[idx]))
            clusters[nearest_idx].append((q, count))
            
        # 2. Update centroids to cluster weighted means
        new_centroids = []
        for idx in range(num_bins):
            cluster = clusters[idx]
            if cluster:
                w_sum = sum(q * count for q, count in cluster)
                w_count = sum(count for q, count in cluster)
                new_centroids.append(w_sum / float(w_count))
            else:
                new_centroids.append(centroids[idx])
                
        # Check convergence
        if all(abs(c1 - c2) < 0.01 for c1, c2 in zip(centroids, new_centroids)):
            break
        centroids = new_centroids
        
    return sorted([int(round(c)) for c in centroids])

class Quantizer:
    """Context-aware Phred quality score quantizer supporting levels 1-5."""
    
    def __init__(self, level: BinningLevel = BinningLevel.LEVEL_3_ADAPTIVE_CONTEXT, lossless: bool = False):
        self.level = level
        self.lossless = lossless or (level == BinningLevel.LEVEL_5_LOSSLESS)

    @staticmethod
    def quantize_illumina8(q: int) -> int:
        if q <= 1: return 0
        if q <= 9: return 6
        if q <= 19: return 15
        if q <= 24: return 22
        if q <= 29: return 27
        if q <= 34: return 33
        if q <= 39: return 37
        return 40

    @staticmethod
    def quantize_coarse4(q: int) -> int:
        if q <= 9: return 6
        if q <= 24: return 18
        if q <= 34: return 30
        return 38

    @staticmethod
    def quantize_binary2(q: int) -> int:
        return 10 if q < 20 else 35

    @classmethod
    def quantize_adaptive(cls, q: int, entropy: float, cycle: float) -> int:
        # Low complexity / homopolymer: protect high fidelity
        if entropy < 1.0:
            return cls.quantize_illumina8(q)
        # High complexity / informative: coarse binning eliminates noise
        if entropy >= 1.5:
            if cycle < 0.75 and 28 <= q <= 34:
                return 30
            return cls.quantize_coarse4(q)
        # Intermediate 6-bin:
        if q <= 2: return 0
        if q <= 12: return 8
        if q <= 22: return 18
        if q <= 29: return 26
        if q <= 36: return 33
        return 39

    def quantize(
        self,
        sequence: str,
        quality_str: str,
        is_read_2: bool = False,
        custom_centroids: List[int] = None
    ) -> Tuple[str, List[int]]:
        """
        Quantizes Phred quality score string given sequence.
        Supports Paired-End Read 1/2 decay asymmetry and dynamic Lloyd-Max centroids.
        Returns (quantized_quality_str, residual_list).
        """
        seq_len = min(len(sequence), len(quality_str))
        quant_chars = []
        residuals = []
        
        for i in range(seq_len):
            raw_q = ord(quality_str[i]) - 33
            if raw_q < 0: raw_q = 0
            
            if custom_centroids:
                # Dynamic Lloyd-Max centroid mapping
                q_hat = min(custom_centroids, key=lambda c: abs(raw_q - c))
            elif self.level == BinningLevel.LEVEL_1_ILLUMINA8:
                q_hat = self.quantize_illumina8(raw_q)
            elif self.level == BinningLevel.LEVEL_2_COARSE4:
                q_hat = self.quantize_coarse4(raw_q)
            elif self.level == BinningLevel.LEVEL_3_ADAPTIVE_CONTEXT or self.level == BinningLevel.LEVEL_5_LOSSLESS:
                ent = calculate_local_entropy(sequence, i, 3)
                cyc = calculate_cycle_decay(i, seq_len, is_read_2=is_read_2)
                q_hat = self.quantize_adaptive(raw_q, ent, cyc)
            elif self.level == BinningLevel.LEVEL_4_BINARY:
                q_hat = self.quantize_binary2(raw_q)
            else:
                q_hat = self.quantize_illumina8(raw_q)
                
            q_hat = min(93, max(0, q_hat))
            quant_chars.append(chr(q_hat + 33))
            
            if self.lossless:
                residuals.append(raw_q - q_hat)
                
        return ''.join(quant_chars), residuals

    @staticmethod
    def restore(quantized_str: str, residuals: List[int]) -> str:
        """Restores exact original Phred quality string using residual stream."""
        if not residuals:
            return quantized_str
            
        n = min(len(quantized_str), len(residuals))
        restored = []
        for i in range(n):
            q_hat = ord(quantized_str[i]) - 33
            q_orig = q_hat + residuals[i]
            q_orig = max(0, min(93, q_orig))
            restored.append(chr(q_orig + 33))
        return ''.join(restored)
