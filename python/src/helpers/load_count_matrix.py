from config import get_config
import boto3
import datetime
import anndata
import os
import tempfile
from helpers.dynamo import get_item_from_dynamo

config = get_config()


def _download_obj(bucket, key):

    tmpdir = tempfile.TemporaryDirectory()

    if config.CLUSTER_ENV == "development" or config.CLUSTER_ENV == "test":
        print("In development, creating a temporary directory...")
        path = tmpdir.name
    else:
        path = "/data"

    try:
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        file_path = os.path.join(path, f"{bucket}_{key}")

        with open(file_path, "wb") as f:
            client.download_fileobj(Bucket=bucket, Key=key, Fileobj=f)

        adata = anndata.read_h5ad(file_path)
        print(datetime.datetime.utcnow(), "File read.")
        return adata

    except Exception as e:
        print(datetime.datetime.utcnow(), "Could not get file from S3", e)
        raise e


def _load_file(matrix_path):
    # intercept here if task is to prepare experiment, don't download from s3
    print(datetime.datetime.utcnow(), "Have to download anndata file from s3")
    bucket, key = matrix_path.split("/", 1)
    adata = _download_obj(bucket, key)
    print(datetime.datetime.utcnow(), "File was loaded.")

    if "cell_ids" not in adata.obs:
        raise ValueError(
            "You must have `cell_ids` in your anndata file for integer cell IDs."
        )

    return adata


def get_adata(adata, experiment_id):
    print(datetime.datetime.utcnow(), "adata does not exist, I need to download it ...")
    matrix_path = get_item_from_dynamo(experiment_id, "matrixPath")
    adata = _load_file(matrix_path)
    return adata
