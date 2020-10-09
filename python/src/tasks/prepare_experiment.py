import scanpy as sc
import boto3
import io
import os
import shutil
import pickle
import json
from result import Result
from helpers.dynamo import get_item_from_dynamo
from config import get_config

config = get_config()


class PrepareExperiment:
    def __init__(self, msg, adata):
        self.directory_path = msg["body"]["sourceMatrixPath"]
        self.source_bucket = msg["body"]["sourceBucket"]
        self.experiment_id = msg["experimentId"]
        self.adata = adata

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps({"url": result})

        # Return a list of formatted results.
        return [Result(result)]

    def _download_file_to_dir(self):
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        objList = client.list_objects(
            Bucket=self.source_bucket, Prefix=self.directory_path
        )
        objList = objList.get("Contents", [])
        if len(objList) == 0:
            raise Exception(
                "Couldn't download count matrix files: path {} in bucket {} is does not exist or is empty".format(
                    self.directory_path, self.source_bucket
                )
            )
        if not os.path.exists(self.directory_path):
            os.makedirs(self.directory_path)
            for obj in objList:
                if obj["Size"] > 0:
                    with open(obj["Key"], "wb") as data:
                        client.download_fileobj(self.source_bucket, obj["Key"], data)

    def _upload_anndata(self, adata):
        matrix_path = get_item_from_dynamo(self.experiment_id, "matrixPath")
        bucket, key = matrix_path.split("/", 1)

        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

        result = io.BytesIO()
        pickle.dump(adata, result)
        result.seek(0)
        client.upload_fileobj(result, bucket, key)
        print("Uploaded Anndata file {} to bucket {}".format(key, bucket))

    def _clean(self):
        if os.path.exists(self.directory_path):
            shutil.rmtree(self.directory_path)
            dir_root = self.directory_path.split("/")[0]
            if os.path.isdir(dir_root):
                os.removedirs(dir_root)
            print("Directory cleaned.")

    def compute(self):
        try:
            self._download_file_to_dir()
            self.adata = sc.read_10x_mtx(path=self.directory_path)
            print("Uploading the new anndata file to s3 ...")
            self._upload_anndata(self.adata)
            self._clean()
            return self._format_result(
                "https://my-experiment-url/{}".format(self.experiment_id)
            )
        except Exception as e:
            self._clean()
            raise e
