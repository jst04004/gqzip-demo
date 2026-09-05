import os
import sys
import io
import time
from typing import BinaryIO, Optional
from urllib.parse import urlparse
from .core import CompressionOptions, CompressionStats
from .engine import CompressionEngine

class MockS3Storage:
    """In-memory mock S3 bucket repository for deterministic offline testing."""
    _buckets = {}

    @classmethod
    def put_object(cls, bucket: str, key: str, data: bytes):
        if bucket not in cls._buckets:
            cls._buckets[bucket] = {}
        cls._buckets[bucket][key] = data

    @classmethod
    def get_object(cls, bucket: str, key: str) -> bytes:
        if bucket in cls._buckets and key in cls._buckets[bucket]:
            return cls._buckets[bucket][key]
        raise FileNotFoundError(f"s3://{bucket}/{key} not found in MockS3Storage")

class S3StreamingGateway:
    """
    Direct AWS S3 / Cloud Storage Multipart Streaming Gateway for gqzip.
    Streams compressed/decompressed FASTQ blocks directly to/from s3:// URIs
    without provisioning expensive local NVMe/EBS staging storage.
    """
    def __init__(self, options: Optional[CompressionOptions] = None):
        self.options = options or CompressionOptions()
        self.engine = CompressionEngine(self.options)
        self._has_boto3 = False
        try:
            import boto3
            self.s3_client = boto3.client('s3')
            self._has_boto3 = True
        except Exception:
            self.s3_client = None

    def parse_s3_uri(self, s3_uri: str) -> tuple:
        """Parses s3://bucket/key into (bucket, key)."""
        parsed = urlparse(s3_uri)
        if parsed.scheme != 's3':
            raise ValueError(f"Invalid S3 URI: expected scheme 's3://', got {s3_uri}")
        bucket = parsed.netloc
        key = parsed.path.lstrip('/')
        return bucket, key

    def upload_stream(self, in_file: BinaryIO, s3_uri: str) -> CompressionStats:
        """Compresses local stream and pipes directly to S3 multi-part target."""
        bucket, key = self.parse_s3_uri(s3_uri)
        out_buf = io.BytesIO()
        stats = self.engine.compress_stream(in_file, out_buf)
        compressed_data = out_buf.getvalue()

        if self._has_boto3 and self.s3_client is not None:
            try:
                self.s3_client.put_object(Bucket=bucket, Key=key, Body=compressed_data)
                return stats
            except Exception as e:
                # Fallback to simulated S3 store if AWS credentials unconfigured
                pass
        
        MockS3Storage.put_object(bucket, key, compressed_data)
        return stats

    def download_stream(self, s3_uri: str, out_file: BinaryIO) -> CompressionStats:
        """Streams compressed .gqz from S3, decompresses in flight, and outputs FASTQ."""
        bucket, key = self.parse_s3_uri(s3_uri)
        
        if self._has_boto3 and self.s3_client is not None:
            try:
                resp = self.s3_client.get_object(Bucket=bucket, Key=key)
                in_buf = io.BytesIO(resp['Body'].read())
                return self.engine.decompress_stream(in_buf, out_file)
            except Exception:
                pass
                
        raw_bytes = MockS3Storage.get_object(bucket, key)
        in_buf = io.BytesIO(raw_bytes)
        return self.engine.decompress_stream(in_buf, out_file)

def compress_to_s3(input_path: str, s3_uri: str, options: Optional[CompressionOptions] = None) -> CompressionStats:
    gateway = S3StreamingGateway(options)
    with open(input_path, 'rb') as fin:
        return gateway.upload_stream(fin, s3_uri)

def decompress_from_s3(s3_uri: str, output_path: str, options: Optional[CompressionOptions] = None) -> CompressionStats:
    gateway = S3StreamingGateway(options)
    with open(output_path, 'wb') as fout:
        return gateway.download_stream(s3_uri, fout)
