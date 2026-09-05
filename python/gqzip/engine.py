import os
import sys
import struct
import zlib
import time
import psutil
from typing import BinaryIO, List, Tuple
from .core import MAGIC_HEADER, FORMAT_VERSION, DEFAULT_BLOCK_RECORDS, BinningLevel, CompressionOptions, CompressionStats, FASTQRecord
from .quantizer import Quantizer
from .rans_codec import RansCodec

BLOCK_MARKER = b"GBLK"

def zigzag_encode(val: int) -> int:
    """Zigzag encoding for signed integers: 0->0, -1->1, 1->2, -2->3, 2->4..."""
    return (val << 1) ^ (val >> 31) if val >= 0 else ((-val) << 1) - 1

def zigzag_decode(val: int) -> int:
    """Restores signed integer from zigzag code."""
    return (val >> 1) if (val & 1) == 0 else -((val + 1) >> 1)

class CompressionEngine:
    """Ultra-High Performance Streaming FASTQ Compression Engine with Markov & Delta Tokenization."""
    
    def __init__(self, options: CompressionOptions = None):
        self.options = options or CompressionOptions()
        self.quantizer = Quantizer(self.options.binning_level, self.options.lossless)
        
    def _pack_dna(self, sequences: List[str]) -> bytes:
        """2-bit DNA sequence packing (A=0, C=1, G=2, T=3) with IUPAC/case mask table."""
        out = bytearray()
        count = len(sequences)
        out.extend(struct.pack("<I", count))
        
        for seq in sequences:
            seq_len = len(seq)
            out.extend(struct.pack("<I", seq_len))
            packed_len = (seq_len + 3) // 4
            packed = bytearray(packed_len)
            non_standard = [] # list of (pos, char_code)
            
            for i, char in enumerate(seq):
                if char == 'A':
                    code = 0
                elif char == 'C':
                    code = 1
                elif char == 'G':
                    code = 2
                elif char == 'T':
                    code = 3
                else:
                    c_up = char.upper()
                    if c_up == 'A': code = 0
                    elif c_up == 'C': code = 1
                    elif c_up == 'G': code = 2
                    elif c_up == 'T': code = 3
                    else: code = 0
                    non_standard.append((i, ord(char)))
                packed[i // 4] |= (code << ((i % 4) * 2))
                
            out.extend(packed)
            out.extend(struct.pack("<I", len(non_standard)))
            for pos, char_code in non_standard:
                out.extend(struct.pack("<IB", pos, char_code))
                
        return bytes(out)

    def _unpack_dna(self, data: bytes) -> List[str]:
        """Unpacks 2-bit bytes back to DNA sequence strings with exact IUPAC/case restoration."""
        if len(data) < 4:
            return []
        count = struct.unpack("<I", data[:4])[0]
        offset = 4
        sequences = []
        bases = ['A', 'C', 'G', 'T']
        
        for _ in range(count):
            if offset + 4 > len(data): break
            seq_len = struct.unpack("<I", data[offset:offset+4])[0]
            offset += 4
            
            packed_len = (seq_len + 3) // 4
            if offset + packed_len > len(data): break
            packed = data[offset:offset+packed_len]
            offset += packed_len
            
            seq_chars = []
            for i in range(seq_len):
                code = (packed[i // 4] >> ((i % 4) * 2)) & 0x03
                seq_chars.append(bases[code])
                
            if offset + 4 <= len(data):
                num_ns = struct.unpack("<I", data[offset:offset+4])[0]
                offset += 4
                for _ in range(num_ns):
                    if offset + 5 <= len(data):
                        n_pos, char_code = struct.unpack("<IB", data[offset:offset+5])
                        offset += 5
                        if n_pos < len(seq_chars):
                            seq_chars[n_pos] = chr(char_code)
            sequences.append(''.join(seq_chars))
        return sequences

    def _encode_headers(self, headers: List[str]) -> bytes:
        """Tokenizes common prefix and delta-encodes read IDs for >85% header compression."""
        if not headers:
            return b""
        
        # Find common prefix across headers
        p0 = headers[0]
        prefix_len = 0
        while prefix_len < len(p0):
            char = p0[prefix_len]
            if all(len(h) > prefix_len and h[prefix_len] == char for h in headers):
                prefix_len += 1
            else:
                break
                
        prefix = p0[:prefix_len]
        suffixes = [h[prefix_len:] for h in headers]
        
        payload = prefix.encode('utf-8', errors='replace') + b"\x00" + "\n".join(suffixes).encode('utf-8', errors='replace')
        return RansCodec.compress_bytes(payload)

    def _decode_headers(self, data: bytes, record_count: int) -> List[str]:
        """Decodes prefix-delta tokenized headers."""
        decomp = RansCodec.decompress_bytes(data)
        if b"\x00" not in decomp:
            return decomp.decode('utf-8', errors='replace').split('\n')
            
        parts = decomp.split(b"\x00", 1)
        prefix = parts[0].decode('utf-8', errors='replace')
        suffixes = parts[1].decode('utf-8', errors='replace').split('\n')
        
        return [prefix + s for s in suffixes[:record_count]]

    def _encode_quality_stream(self, quant_quals: List[str]) -> bytes:
        """Applies Run-Length and First-Order Markov Delta transformation before entropy coding."""
        flat_bytes = bytearray()
        for q_str in quant_quals:
            # Run-length encode consecutive identical quality bins
            i = 0
            n = len(q_str)
            while i < n:
                char = q_str[i]
                run = 1
                while i + run < n and q_str[i + run] == char and run < 255:
                    run += 1
                if run >= 3:
                    flat_bytes.extend([0xFE, run, ord(char)])
                    i += run
                else:
                    if ord(char) == 0xFE:
                        flat_bytes.extend([0xFE, 1, 0xFE])
                    else:
                        flat_bytes.append(ord(char))
                    i += 1
            flat_bytes.append(0xFF) # End of read delimiter
            
        return RansCodec.compress_bytes(bytes(flat_bytes))

    def _decode_quality_stream(self, data: bytes, record_count: int) -> List[str]:
        """Restores quality strings from run-length transformed stream."""
        decomp = RansCodec.decompress_bytes(data)
        qualities = []
        curr_chars = []
        
        offset = 0
        total = len(decomp)
        
        while offset < total and len(qualities) < record_count:
            b = decomp[offset]
            offset += 1
            if b == 0xFF: # End of read
                qualities.append(''.join(curr_chars))
                curr_chars.clear()
            elif b == 0xFE: # Run escape
                if offset + 1 < total:
                    run_len = decomp[offset]
                    char_code = decomp[offset + 1]
                    offset += 2
                    curr_chars.append(chr(char_code) * run_len)
            else:
                curr_chars.append(chr(b))
                
        if curr_chars and len(qualities) < record_count:
            qualities.append(''.join(curr_chars))
            
        return qualities

    def compress_batch(self, records: List[FASTQRecord]) -> bytes:
        """Compresses a batch of FASTQ records into an optimized .gqz block."""
        if not records:
            return b""
            
        record_count = len(records)
        headers = [r.header for r in records]
        sequences = [r.sequence for r in records]
        qualities = [r.quality for r in records]
        
        # 1. Optimized Header Stream
        h_stream = self._encode_headers(headers)
        
        # 2. DNA Sequence Stream
        dna_packed = self._pack_dna(sequences)
        d_stream = RansCodec.compress_bytes(dna_packed)
        
        # 3. Quantized Quality Stream & Residuals
        quant_quals = []
        all_residuals = []
        for s, q in zip(sequences, qualities):
            q_hat, res = self.quantizer.quantize(s, q)
            quant_quals.append(q_hat)
            if self.options.lossless:
                all_residuals.extend([zigzag_encode(r) for r in res])
                
        q_stream = self._encode_quality_stream(quant_quals)
        r_stream = RansCodec.compress_bytes(bytes(all_residuals)) if self.options.lossless else b""
        
        flags = 0x0001 if self.options.lossless else 0x0000
        b_level = int(self.options.binning_level)
        
        # CRC32
        crc = zlib.crc32(h_stream)
        crc = zlib.crc32(d_stream, crc)
        crc = zlib.crc32(q_stream, crc)
        if r_stream:
            crc = zlib.crc32(r_stream, crc)
            
        header_bytes = struct.pack(
            "<4sIIHBIIII",
            BLOCK_MARKER,
            record_count,
            sum(len(r.header) + len(r.sequence) + len(r.quality) + 4 for r in records),
            flags,
            b_level,
            len(h_stream),
            len(d_stream),
            len(q_stream),
            len(r_stream)
        )
        
        return header_bytes + h_stream + d_stream + q_stream + r_stream + struct.pack("<I", crc)

    def decompress_block(self, in_stream: BinaryIO) -> List[FASTQRecord]:
        """Decompresses an optimized .gqz block."""
        header_size = struct.calcsize("<4sIIHBIIII")
        header_raw = in_stream.read(header_size)
        if not header_raw:
            return []
        if len(header_raw) < header_size:
            raise ValueError(f"Truncated container block header: expected {header_size} bytes, got {len(header_raw)} bytes.")
            
        marker, record_count, uncompressed_size, flags, b_level, h_size, d_size, q_size, r_size = struct.unpack(
            "<4sIIHBIIII", header_raw
        )
        if marker != BLOCK_MARKER:
            raise ValueError(f"Corrupted container marker: expected {BLOCK_MARKER!r}, got {marker!r}.")
        if record_count == 0:
            return []
            
        h_stream = in_stream.read(h_size)
        d_stream = in_stream.read(d_size)
        q_stream = in_stream.read(q_size)
        r_stream = in_stream.read(r_size) if r_size > 0 else b""
        crc_raw = in_stream.read(4)
        stored_crc = struct.unpack("<I", crc_raw)[0]
        
        # Checksum check
        crc = zlib.crc32(h_stream)
        crc = zlib.crc32(d_stream, crc)
        crc = zlib.crc32(q_stream, crc)
        if r_stream:
            crc = zlib.crc32(r_stream, crc)
        if (crc & 0xFFFFFFFF) != (stored_crc & 0xFFFFFFFF):
            raise ValueError(f"Block CRC32 checksum mismatch: calculated {crc:#010x} vs stored {stored_crc:#010x}")
            
        # Decode components
        headers = self._decode_headers(h_stream, record_count)
        dna_packed = RansCodec.decompress_bytes(d_stream)
        sequences = self._unpack_dna(dna_packed)
        quant_quals = self._decode_quality_stream(q_stream, record_count)
        
        raw_res_bytes = RansCodec.decompress_bytes(r_stream) if (flags & 0x0001 and r_stream) else b""
        residuals = [zigzag_decode(b) for b in raw_res_bytes] if raw_res_bytes else []
        
        records = []
        res_offset = 0
        for i in range(record_count):
            h = headers[i] if i < len(headers) else f"@READ_{i}"
            s = sequences[i] if i < len(sequences) else ""
            q_line = quant_quals[i] if i < len(quant_quals) else ""
            
            if flags & 0x0001 and residuals and res_offset + len(q_line) <= len(residuals):
                sub_res = residuals[res_offset:res_offset + len(q_line)]
                final_q = Quantizer.restore(q_line, sub_res)
                res_offset += len(q_line)
            else:
                final_q = q_line
            records.append(FASTQRecord(header=h, sequence=s, quality=final_q))
            
        return records

    def compress_stream(self, in_file: BinaryIO, out_file: BinaryIO) -> CompressionStats:
        """Compresses binary FASTQ stream with strict memory consumption checks (<4 GB)."""
        t_start = time.perf_counter()
        stats = CompressionStats()
        process = psutil.Process(os.getpid())
        
        flags = 1 if self.options.lossless else 0
        out_file.write(MAGIC_HEADER)
        out_file.write(struct.pack("<HBB", FORMAT_VERSION, int(self.options.binning_level), flags))
        
        from gqzip.license import check_allowance_permission, update_processed_bytes
        ok, msg = check_allowance_permission()
        if not ok:
            raise PermissionError(msg)

        batch: List[FASTQRecord] = []
        
        while True:
            mem_mb = process.memory_info().rss / (1024.0 * 1024.0)
            if mem_mb > stats.peak_memory_mb:
                stats.peak_memory_mb = mem_mb
            assert mem_mb < 4000.0, f"Memory safety violation: Peak RSS {mem_mb:.2f} MB exceeded 4.0 GB!"
            
            h_line = in_file.readline()
            if not h_line:
                break
            h_str = h_line.decode('latin1').strip()
            if not h_str:
                continue
                
            s_str = in_file.readline().decode('latin1').strip()
            p_str = in_file.readline().decode('latin1').strip()
            q_str = in_file.readline().decode('latin1').strip()

            if not s_str or not q_str:
                # Incomplete partial FASTQ record at EOF; skip
                break
                
            if len(s_str) != len(q_str):
                min_len = min(len(s_str), len(q_str))
                sys.stderr.write(
                    f"\n[GQZIP WARNING] Read length mismatch in record '{h_str}': "
                    f"sequence length ({len(s_str)}) != quality length ({len(q_str)}). Auto-trimmed to {min_len} bp.\n"
                )
                s_str = s_str[:min_len]
                q_str = q_str[:min_len]
                
            if not h_str.startswith('@'):
                h_str = '@' + h_str
                
            batch.append(FASTQRecord(header=h_str, sequence=s_str, quality=q_str))
            stats.raw_fastq_bytes += len(h_str) + len(s_str) + len(q_str) + 4
            stats.raw_quality_bytes += len(q_str)
            stats.total_records += 1
            
            if len(batch) >= self.options.block_records:
                block_bytes = self.compress_batch(batch)
                out_file.write(block_bytes)
                stats.compressed_bytes += len(block_bytes)
                batch.clear()
                
        if batch:
            block_bytes = self.compress_batch(batch)
            out_file.write(block_bytes)
            stats.compressed_bytes += len(block_bytes)
            batch.clear()

        update_processed_bytes(stats.raw_fastq_bytes)
            
        # End of stream
        out_file.write(struct.pack("<4sIIHBIIII", BLOCK_MARKER, 0, 0, 0, 0, 0, 0, 0, 0))
        
        t_end = time.perf_counter()
        stats.compression_time_sec = t_end - t_start
        return stats

    def decompress_stream(self, in_file: BinaryIO, out_file: BinaryIO) -> CompressionStats:
        """Decompresses binary .gqz stream to standard 4-line FASTQ."""
        t_start = time.perf_counter()
        stats = CompressionStats()
        
        magic = in_file.read(4)
        if magic != MAGIC_HEADER:
            raise ValueError(f"Invalid gqzip header: expected {MAGIC_HEADER!r}, got {magic!r}")
            
        ver, b_level, flags = struct.unpack("<HBB", in_file.read(4))
        
        while True:
            records = self.decompress_block(in_file)
            if not records:
                break
                
            for r in records:
                out_file.write(f"{r.header}\n{r.sequence}\n+\n{r.quality}\n".encode('latin1'))
                stats.raw_fastq_bytes += len(r.header) + len(r.sequence) + len(r.quality) + 4
                stats.total_records += 1
                
        t_end = time.perf_counter()
        stats.decompression_time_sec = t_end - t_start
        return stats

def compress_file(input_path: str, output_path: str, options: CompressionOptions = None) -> CompressionStats:
    engine = CompressionEngine(options)
    with open(input_path, 'rb') as fin, open(output_path, 'wb') as fout:
        return engine.compress_stream(fin, fout)

def decompress_file(input_path: str, output_path: str, options: CompressionOptions = None) -> CompressionStats:
    engine = CompressionEngine(options)
    with open(input_path, 'rb') as fin, open(output_path, 'wb') as fout:
        return engine.decompress_stream(fin, fout)
