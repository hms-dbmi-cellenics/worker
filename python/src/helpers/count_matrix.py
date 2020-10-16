from config import get_config
import boto3
import datetime
import os
import hashlib
from helpers.dynamo import get_item_from_dynamo

config = get_config()


def _download_obj(bucket, key):
    try:
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        print("about to download file ", bucket, key)
        experiment_files = client.list_objects(
            Bucket="biomage-source-development", Prefix=key
        )
        path = f"/data/{key}"

        print("MY PATH:    ", path)
        if not os.path.exists(path):
            print("not existing, creating: ")
            os.makedirs(path)
        for content in experiment_files["Contents"]:
            print("RESP     ", content)
            file_name = content["Key"].split("/", 1)[-1]
            print("file name: ", file_name)
            full_path = f"{path}/{file_name}"
            print("FULL PATH:     ", full_path)
            with open(full_path, "wb+") as f:
                client.download_fileobj(
                    Bucket="biomage-source-development",
                    Key=content["Key"],
                    Fileobj=f,
                )
                f.seek(0)
    except Exception as e:
        print(datetime.datetime.utcnow(), "Could not get file from S3", e)
        raise e
    print(datetime.datetime.utcnow(), "File was loaded.")


def _get_file_name(experiment_id):
    file_name = "python.h5"
    path = f"/data/{experiment_id}"
    adata_file = f"{path}/{file_name}"
    return adata_file


def get_adata_path(experiment_id):
    print(
        datetime.datetime.utcnow(),
        "adata does not exist or has changed, I need to download it ...",
    )
    adata_path = _get_file_name(experiment_id)
    _download_obj(config.BUCKET_NAME, experiment_id)
    return adata_path


def calculate_file_etag(file_path, chunk_size=8 * 1024 * 1024):
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


def is_file_changed(adata, experiment_id, adata_path):
    client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
    key = adata_path.strip("/data")
    print("key:    ", key)
    resp = client.head_object(Bucket="biomage-source-development", Key=key)
    s3_etag = resp["ETag"]
    file_etag = calculate_file_etag(adata_path)

    return s3_etag != file_etag
