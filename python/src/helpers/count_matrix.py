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

        objects = {o["Key"]: o["ETag"] for o in objects if o["Size"] > 0}

        return objects

    def md5_full_file(self, file_path):
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def megabyte_chunks(self, filesize, num_parts):
        x = filesize / int(num_parts)
        y = x % 1048576
        return int(x + 1048576 - y)

    def validate_etag(self, file_path, etag):

        # get the size of the file. if it doesn't exist,
        # the etags definitely don't match
        try:
            filesize = os.path.getsize(file_path)
        except OSError:
            return False

        # get the number of parts. if this doesn't exist,
        # the etag has no parts, so it is just the md5 digest
        etag = etag.replace('"', "")

        try:
            num_parts = int(etag.split("-")[1])
        except IndexError:
            return etag == self.md5_full_file(file_path)

        # these are some common part sizes used. the last one
        # is the file size divided by the number of parts,
        # rounded to the nearest megabyte.
        partsizes = [
            8388608,
            15728640,
            self.megabyte_chunks(filesize, num_parts),
        ]

        valid_sizes = lambda partsize: (
            partsize < filesize and (float(filesize) / float(partsize)) <= num_parts
        )

        for partsize in filter(valid_sizes, partsizes):
            md5_digests = []

            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(partsize), b""):
                    md5_digests.append(hashlib.md5(chunk).digest())

            candidate = (
                hashlib.md5(b"".join(md5_digests)).hexdigest()
                + "-"
                + str(len(md5_digests))
            )

            if candidate == etag:
                return True

        return False

    def download_object(self, key, etag):
        path = os.path.join(config.LOCAL_DIR, key)

        print(f"Now checking {key} (etag: {etag}) in S3 against {path}...")

        if self.path_exists and self.validate_etag(path, etag):
            print("Skipping downloads as etags match.")
            return False

        print(f"Downloading {key} (etag: {etag})...")

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

        print(
            datetime.datetime.utcnow(),
            "Found",
            len(objects),
            "objects matching experiment.",
        )

        synced = {key: self.download_object(key, etag) for key, etag in objects.items()}
        adata_path = os.path.join("../../data/test/Objeto", "python.h5ad")
        # adata_path = os.path.join(self.local_path, "data/test/Objeto", "python.h5ad")
        self.update_anndata(synced, adata_path)
