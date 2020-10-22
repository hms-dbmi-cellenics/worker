from config import get_config
import boto3
import datetime
import os
import hashlib
import anndata
from config import get_config

config = get_config()


class CountMatrix:
    def __init__(self):
        self.config = get_config()
        self.local_path = os.path.join(self.config.LOCAL_DIR, self.config.EXPERIMENT_ID)
        self.s3 = boto3.client("s3", **self.config.BOTO_RESOURCE_KWARGS)

        self.adata = None

    def get_objects(self):
        objects = self.s3.list_objects_v2(
            Bucket=self.config.SOURCE_BUCKET, Prefix=self.config.EXPERIMENT_ID
        )
        objects = objects.get("Contents")
        if not objects:
            return {}
        objects = {o["Key"]: o["ETag"] for o in objects}

        return objects

    def calculate_file_etag(self, file_path, chunk_size=8 * 1024 * 1024):
        md5s = []

        try:
            with open(file_path, "rb") as fp:
                while True:
                    data = fp.read(chunk_size)
                    if not data:
                        break
                    md5s.append(hashlib.md5(data))
        except OSError as e:
            pass

        if len(md5s) < 1:
            return '"{}"'.format(hashlib.md5().hexdigest())

        if len(md5s) == 1:
            return '"{}"'.format(md5s[0].hexdigest())

        digests = b"".join(m.digest() for m in md5s)
        digests_md5 = hashlib.md5(digests)
        return '"{}-{}"'.format(digests_md5.hexdigest(), len(md5s))

    def download_object(self, key, etag):
        path = os.path.join(config.LOCAL_DIR, key)

        if self.path_exists:
            print("We have files from previous runs, comparing checksums to etag...")
            local_etag = self.calculate_file_etag(path)

            if local_etag == etag:
                print("Skipping downloads as etags match.")
                return False

        print(f"Downloading {key} (etag: {etag}) from S3 to {path}...")

        with open(path, "wb+") as f:
            self.s3.download_fileobj(
                Bucket=self.config.SOURCE_BUCKET,
                Key=key,
                Fileobj=f,
            )

            f.seek(0)

        return True

    def update_anndata(self, synced, adata_path):
        adata_key = os.path.join(self.config.EXPERIMENT_ID, "python.h5ad")

        if not self.adata:
            print("AnnData does not exist in memory, creating it now...")
            self.adata = anndata.read_h5ad(adata_path)
        elif self.adata and synced[adata_key]:
            print("AnnData has been updated, reloading...")
            self.adata = anndata.read_h5ad(adata_path)
        else:
            raise ValueError(
                "AnnData could not be loaded due to an unknown error."
            )

        if "cell_ids" not in self.adata.obs:
            raise ValueError(
                "You must have `cell_ids` in your anndata file for integer cell IDs."
            )

    def sync(self):
        # check if path existed before running this
        self.path_exists = os.path.exists(self.local_path)

        if not self.path_exists:
            print(
                datetime.datetime.utcnow(),
                "Path",
                self.local_path,
                "does not yet exist, creating it...",
            )
            os.makedirs(self.local_path)

        # get object in bucket and their etags
        objects = self.get_objects()

        synced = {key: self.download_object(key, etag) for key, etag in objects.items()}

        adata_path = os.path.join(self.local_path, "python.h5ad")
        self.update_anndata(synced, adata_path)
