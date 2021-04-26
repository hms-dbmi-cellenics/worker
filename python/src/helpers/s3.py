import boto3
import json
from config import get_config
import os

import aws_xray_sdk as xray

config = get_config()


def get_cell_sets(experiment_id):
    dir_path = os.path.join(config.LOCAL_DIR, f"{experiment_id}")

    path_exists = os.path.exists(dir_path)

    if not path_exists:
        os.makedirs(dir_path)

    with open(f"{dir_path}/cell_sets.json", "wb+") as f:
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        print("downloading cellsets")

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