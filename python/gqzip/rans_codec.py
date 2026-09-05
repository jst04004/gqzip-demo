"""
Pure 32-bit Asymmetric Numeral Systems (rANS) Entropy Codec Layer with Order-1 Markov Context Modeling.
Supports high-performance byte-stream and residual entropy coding. Zero zstd dependencies.
"""

import struct
from typing import List

# Try importing native C++ SIMD rANS extension module
try:
    import _gqzip_cpp
    _HAS_CPP = True
except ImportError:
    _HAS_CPP = False

RANS_SCALE_BITS = 12
RANS_TOTAL_FREQ = 1 << RANS_SCALE_BITS  # 4096
RANS_BYTE_L = 1 << 23  # Normalized lower bound

class RansCodec:
    """Interleaved 32-bit rANS Entropy Codec with Order-1 Markov Context Modeling."""

    @classmethod
    def compress_bytes(cls, data: bytes) -> bytes:
        if not data:
            return b""
        if _HAS_CPP:
            return _gqzip_cpp.RansCodec.encode_order1(data)
        return cls._encode_python_rans(data)

    @classmethod
    def decompress_bytes(cls, compressed: bytes, max_output_size: int = 100 * 1024 * 1024) -> bytes:
        if not compressed:
            return b""
        if _HAS_CPP:
            return _gqzip_cpp.RansCodec.decode_order1(compressed, max_output_size)
        return cls._decode_python_rans(compressed)

    @classmethod
    def compress_residuals(cls, residuals: List[int]) -> bytes:
        if not residuals:
            return b""
        zigzag = bytearray()
        for r in residuals:
            val = (r << 1) ^ (r >> 31) if r >= 0 else ((-r) << 1) - 1
            zigzag.append(val & 0xFF)
        return cls.compress_bytes(bytes(zigzag))

    @classmethod
    def decompress_residuals(cls, compressed: bytes) -> List[int]:
        if not compressed:
            return []
        raw_bytes = cls.decompress_bytes(compressed)
        res = []
        for b in raw_bytes:
            val = (b >> 1) if (b & 1) == 0 else -((b + 1) >> 1)
            res.append(val)
        return res

    @classmethod
    def _encode_python_rans(cls, data: bytes) -> bytes:
        """Pure Python 32-bit rANS streaming encoder with 12-bit frequency table header."""
        n = len(data)
        if n == 0:
            return b""
            
        counts = [0] * 256
        for b in data:
            counts[b] += 1

        active = [i for i in range(256) if counts[i] > 0]
        if not active:
            return b""

        freqs = [0] * 256
        rem = RANS_TOTAL_FREQ - len(active)
        for i in active:
            freqs[i] = 1 + (counts[i] * rem) // n

        freqs[active[0]] += RANS_TOTAL_FREQ - sum(freqs)

        cum = [0] * 257
        for i in range(256):
            cum[i + 1] = cum[i] + freqs[i]

        x = RANS_BYTE_L
        buf = bytearray()

        for b in reversed(data):
            f = freqs[b]
            c = cum[b]
            x_max = ((RANS_BYTE_L >> 12) << 8) * f
            while x >= x_max:
                buf.append(x & 0xFF)
                x >>= 8
            x = ((x // f) << 12) + (x % f) + c

        payload = struct.pack("<I", x) + bytes(reversed(buf))
        header = struct.pack("<I", n) + struct.pack("<256H", *freqs)
        return header + payload

    @classmethod
    def _decode_python_rans(cls, compressed: bytes) -> bytes:
        """Pure Python 32-bit rANS streaming decoder."""
        header_size = 4 + 512
        if len(compressed) < header_size + 4:
            return b""

        orig_size = struct.unpack("<I", compressed[:4])[0]
        freqs = list(struct.unpack("<256H", compressed[4:header_size]))

        cum = [0] * 257
        for i in range(256):
            cum[i + 1] = cum[i] + freqs[i]

        lookup = [0] * RANS_TOTAL_FREQ
        for i in range(256):
            if freqs[i] > 0:
                for slot in range(cum[i], cum[i + 1]):
                    lookup[slot] = i

        payload = compressed[header_size:]
        if len(payload) < 4:
            return b""

        x_dec = struct.unpack("<I", payload[:4])[0]
        ptr = 4
        output = bytearray()

        for _ in range(orig_size):
            slot = x_dec & (RANS_TOTAL_FREQ - 1)
            sym = lookup[slot]
            output.append(sym)

            f = freqs[sym]
            c = cum[sym]
            x_dec = f * (x_dec >> 12) + (slot - c)

            while x_dec < RANS_BYTE_L and ptr < len(payload):
                x_dec = (x_dec << 8) | payload[ptr]
                ptr += 1

        return bytes(output)
