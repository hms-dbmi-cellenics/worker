from config import get_config
import boto3
import datetime
import os
import hashlib
from helpers.dynamo import get_item_from_dynamo

config = get_config()


# def _download_obj(bucket, key, experiment_id):
#     try:
#         client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
#         print("about to download file ")
#         with tempfile.TemporaryFile(mode="w+b") as f:
#             client.download_fileobj(Bucket=bucket, Key=key, Fileobj=f)
#             f.seek(0)
#             # adata = anndata.read_h5ad(f)
#     except Exception as e:
#         print(datetime.datetime.utcnow(), "Could not get file from S3", e)
#         raise e
#     print(datetime.datetime.utcnow(), "File was loaded.")
#     return adata

# def _save_file_to_disk(adata, experiment_id, key):
#     file_name = key.split(".")[0]

#     path = f"/data/{experiment_id}"
#     print("MY PATH: ", path)
#     if not os.path.exists(path):
#         os.makedirs(path)
#     adata_file = f"{path}/{file_name}.h5ad"
#     print("MY FILE: ", adata_file)
#     adata.write(filename=adata_file)
#     print("Adata file written successfully to disk.")
#     return path


def _download_obj(bucket, key, filename):
    try:
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        print("about to download file ")
        with open(filename, "wb+") as f:
            client.download_fileobj(Bucket=bucket, Key=key, Fileobj=f)
            f.seek(0)
    except Exception as e:
        print(datetime.datetime.utcnow(), "Could not get file from S3", e)
        raise e
    print(datetime.datetime.utcnow(), "File was loaded.")


def _get_file_name(experiment_id, key):
    file_name = key.split(".")[0]
    path = f"/data/{experiment_id}"
    if not os.path.exists(path):
        os.makedirs(path)
    adata_file = f"{path}/{file_name}.h5ad"
    return adata_file


def get_adata_path(experiment_id):
    print(
        datetime.datetime.utcnow(),
        "adata does not exist or has changed, I need to download it ...",
    )
    matrix_path = get_item_from_dynamo(experiment_id, "matrixPath")
    bucket, key = matrix_path.split("/", 1)
    adata_path = _get_file_name(experiment_id, key)
    _download_obj(bucket, key, adata_path)
    return adata_path


def is_file_changed(adata, experiment_id):
    # compare hashes
    matrix_path = get_item_from_dynamo(experiment_id, "matrixPath")
    bucket, key = matrix_path.split("/", 1)
    client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
    resp = client.head_object(Bucket=bucket, Key=key)
    etag = resp["ETag"].strip('"')

    file_etag = hashlib.md5(open("/iva.h5ad", "rb").read()).hexdigest()
    print("ETAG: ", etag)
    print("file etag: ", file_etag)
    return False
