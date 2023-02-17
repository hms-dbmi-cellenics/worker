import gzip
import json
import os
from logging import info
from pathlib import Path

import aws_xray_sdk as xray
import boto3

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
        cell_sets = json.loads(cell_sets_string)

        return cell_sets["cellSets"]

def get_cell_sets_dict(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        if not cell_sets_dict.get(cell_class["key"]):
          cell_sets_dict[cell_class["key"]] = {
              "key": cell_class["key"],
              "name": cell_class["name"],
              "rootNode": True,
              "childrenKeys": [],
          }

        for cell_set in cell_class["children"]:
            cell_sets_dict[cell_class["key"]]["childrenKeys"].append(cell_set["key"])
            cell_sets_dict[cell_set["key"]] = cell_set

    return cell_sets_dict


# Get all cell sets that match the subset_keys
def get_cell_ids(subset_keys, cell_sets_dict):
    cell_ids = set()

    for subset_key in subset_keys:

        # If cell class, get cell ids from child keys
        if cell_sets_dict[subset_key]["rootNode"]:
            cell_ids = cell_ids.union(
                get_cell_ids(cell_sets_dict[subset_key]["childrenKeys"], cell_sets_dict)
            )
        else:
          cell_ids = cell_ids.union(set(cell_sets_dict[subset_key]['cellIds']))

    return cell_ids


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

        embedding_string = gzip.decompress(f.read())
        embedding = json.loads(embedding_string)

        if(format_for_r):
          # NULL values are deleted in R objects whereas NAs are an indicator of a missing value
          embedding = [ e if e is not None else ["NA", "NA"] for e in embedding ]

        return embedding
