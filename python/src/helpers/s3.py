import boto3
import json
from config import get_config
import os

from helpers.dynamo import get_item_from_dynamo

config = get_config()


def get_cell_sets(experiment_id):
    dir_path = os.path.join(config.LOCAL_DIR, f"{experiment_id}")

    path_exists = os.path.exists(dir_path)

    if not path_exists:
        os.makedirs(dir_path)

    with open(f"{dir_path}/cell_sets.json", "wb+") as f:
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        print("downloading cellsets")

        s3.download_fileobj(
            Bucket=config.CELL_SETS_BUCKET,
            Key=experiment_id,
            Fileobj=f,
        )

        f.seek(0)

        cell_sets_string = f.read()

        cell_sets = json.loads(cell_sets_string)

        return cell_sets["cellSets"]