import os
import json
from logging import info
from pathlib import Path

import boto3
from config import config

import aws_xray_sdk as xray




def get_cell_sets(experiment_id):
    dir_path = os.path.join(config.LOCAL_DIR, f"{experiment_id}")

    Path(dir_path).mkdir(parents=True, exist_ok=True)

    with open(f"{dir_path}/cell_sets.json", "wb+") as f:
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        info("downloading cellsets")

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