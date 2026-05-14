import gzip
import orjson
import os
from logging import info
from pathlib import Path

import aws_xray_sdk as xray
import boto3
import zstandard as zstd

from ..config import config


def get_cell_sets(experiment_id):
    dir_path = os.path.join(config.LOCAL_DIR, f"{experiment_id}")

    Path(dir_path).mkdir(parents=True, exist_ok=True)

    with open(f"{dir_path}/cell_sets.json", "wb+") as f:
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        info(f"Downloading cellsets for experiment {experiment_id}")

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        s3.download_fileobj(
            Bucket=config.CELL_SETS_BUCKET,
            Key=experiment_id,
            Fileobj=f,
        )

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        f.seek(0)

        cell_sets_string = f.read()
        cell_sets = orjson.loads(cell_sets_string)

        return cell_sets["cellSets"]


def get_embedding(etag, format_for_r):
    dir_path = os.path.join(config.LOCAL_DIR, f"{etag}")

    Path(dir_path).mkdir(parents=True, exist_ok=True)

    with open(f"{dir_path}/embedding.json", "wb+") as f:
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        info(f"Downloading embedding with ETag {etag}")

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        s3.download_fileobj(
            Bucket=config.RESULTS_BUCKET,
            Key=etag,
            Fileobj=f,
        )

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        f.seek(0)

        # Try zstd decompression first, fall back to gzip for backwards compatibility
        compressed_data = f.read()
        try:
            dctx = zstd.ZstdDecompressor()
            embedding_string = dctx.decompress(compressed_data)
            info("Decompressed embedding with zstd")
        except zstd.ZstdError:
            info("Failed to decompress with zstd, trying gzip")
            embedding_string = gzip.decompress(compressed_data)
        
        embedding = orjson.loads(embedding_string)

        if(format_for_r):
          # NULL values are deleted in R objects whereas NAs are an indicator of a missing value
          embedding = [ e if e is not None else ["NA", "NA"] for e in embedding ]

        return embedding
