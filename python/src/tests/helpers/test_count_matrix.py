import boto3
from moto import mock_s3

import shutil
from pathlib import Path
from helpers.count_matrix import CountMatrix
from config import get_config

config = get_config()


class TestCountMatrix:
    def get_count_matrix_instance(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        bucket = config.SOURCE_BUCKET
        self.key = f"{config.EXPERIMENT_ID}/python.h5ad"
        s3.create_bucket(
            Bucket=bucket,
            CreateBucketConfiguration={"LocationConstraint": config.AWS_REGION},
        )

        with open("tests/test.h5ad", "rb") as f:
            s3.upload_fileobj(Fileobj=f, Bucket=bucket, Key=self.key)

        self.count_matrix = CountMatrix()
        self.count_matrix.s3 = s3

    @mock_s3
    def test_get_objects(self):
        self.get_count_matrix_instance()
        objs = self.count_matrix.get_objects()
        assert objs == {
            "5e959f9c9f4b120771249001/python.h5ad": '"ad429334fb2de1b0eb8077ba7e222941-5"'
        }

    @mock_s3
    def test_download_object_does_not_exist(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = False
        Path(f"/data/{config.EXPERIMENT_ID}").mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(self.key, "567")
        assert is_downloaded

    @mock_s3
    def test_download_object_previously_existing_etag_doesnt_match(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        Path(f"/data/{config.EXPERIMENT_ID}").mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(self.key, "567")
        assert is_downloaded

    @mock_s3
    def test_download_object_previously_existing_etag_does_match(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        Path(f"/data/{config.EXPERIMENT_ID}").mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(
            self.key, self.count_matrix.calculate_file_etag("tests/test.h5ad")
        )
        assert not is_downloaded

    @mock_s3
    def test_sync_no_previous_data(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        shutil.rmtree(f"/data/{config.EXPERIMENT_ID}", ignore_errors=True)

        self.count_matrix.sync()
        assert "AnnData" in type(self.count_matrix.adata).__name__
        assert not self.count_matrix.path_exists
        assert Path(f"/data/{config.EXPERIMENT_ID}").exists()
        assert self.count_matrix.calculate_file_etag(
            f"/data/{self.key}"
        ) == self.count_matrix.calculate_file_etag("tests/test.h5ad")

    @mock_s3
    def test_sync_previous_data_changed(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True

        shutil.rmtree(f"/data/{config.EXPERIMENT_ID}", ignore_errors=True)
        Path(f"/data/{config.EXPERIMENT_ID}").mkdir(parents=True, exist_ok=True)
        Path(f"/data/{self.key}").touch(exist_ok=True)

        assert self.count_matrix.calculate_file_etag(
            f"/data/{self.key}"
        ) != self.count_matrix.calculate_file_etag("tests/test.h5ad")

        self.count_matrix.sync()
        assert "AnnData" in type(self.count_matrix.adata).__name__
        assert self.count_matrix.path_exists
        assert Path(f"/data/{config.EXPERIMENT_ID}").exists()
        assert self.count_matrix.calculate_file_etag(
            f"/data/{self.key}"
        ) == self.count_matrix.calculate_file_etag("tests/test.h5ad")

    @mock_s3
    def test_sync_previous_data_not_changed(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True

        Path(f"/data/{config.EXPERIMENT_ID}").mkdir(parents=True, exist_ok=True)
        shutil.copy("tests/test.h5ad", f"/data/{self.key}")

        self.count_matrix.sync()
        assert "AnnData" in type(self.count_matrix.adata).__name__
        assert self.count_matrix.path_exists
        assert Path(f"/data/{config.EXPERIMENT_ID}").exists()
        assert self.count_matrix.calculate_file_etag(
            f"/data/{self.key}"
        ) == self.count_matrix.calculate_file_etag("tests/test.h5ad")
