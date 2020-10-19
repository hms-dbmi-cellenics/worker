from config import get_config
import boto3
import datetime
import os
import hashlib

config = get_config()

PATH_TO_FILES = f"/data/{config.EXPERIMENT_ID}"
ADATA_FILE_NAME = "python.h5ad"


def _create_path_to_files():
    if not os.path.exists(PATH_TO_FILES):
        print("not existing, creating: ")
        os.makedirs(PATH_TO_FILES)


def _calculate_file_etag(file_path, chunk_size=8 * 1024 * 1024):
    md5s = []

    with open(file_path, "rb") as fp:
        while True:
            data = fp.read(chunk_size)
            if not data:
                break
            md5s.append(hashlib.md5(data))

    if len(md5s) < 1:
        return '"{}"'.format(hashlib.md5().hexdigest())

    if len(md5s) == 1:
        return '"{}"'.format(md5s[0].hexdigest())

    digests = b"".join(m.digest() for m in md5s)
    digests_md5 = hashlib.md5(digests)
    return '"{}-{}"'.format(digests_md5.hexdigest(), len(md5s))


def is_file_changed(adata, adata_path):
    print("FILE NAME: ", adata_path)
    client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
    file_name = adata_path.split("/")[-1]

    print("key:    ", file_name)
    resp = client.head_object(
        Bucket=config.BUCKET_NAME, Key=f"{config.EXPERIMENT_ID}/{file_name}"
    )
    s3_etag = resp["ETag"]
    file_etag = _calculate_file_etag(adata_path)

    return s3_etag != file_etag


def get_adata_path():
    return f"{PATH_TO_FILES}/{ADATA_FILE_NAME}"


def download_all_files():
    print(
        datetime.datetime.utcnow(),
        "Downloading all files from S3 ...",
    )
    _create_path_to_files()
    try:
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        all_files = client.list_objects(
            Bucket=config.BUCKET_NAME, Prefix=config.EXPERIMENT_ID
        )
        for file_info in all_files["Contents"]:
            file_name = file_info["Key"].split("/", 1)[-1]
            full_path = f"{PATH_TO_FILES}/{file_name}"
            print(f"Downloading file {file_name} to {full_path}")
            with open(full_path, "wb+") as f:
                client.download_fileobj(
                    Bucket=config.BUCKET_NAME,
                    Key=file_info["Key"],
                    Fileobj=f,
                )
                f.seek(0)
    except Exception as e:
        print(datetime.datetime.utcnow(), "Could not get files from S3: ", e)
        raise e
    print(datetime.datetime.utcnow(), "Files were downloaded.")