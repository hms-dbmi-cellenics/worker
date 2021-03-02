import os
import shutil
import boto3
from moto import mock_s3
from pathlib import Path
from helpers.count_matrix import CountMatrix
from config import get_config
import hashlib

config = get_config()

def md5_for_file(path, block_size=256*128, hr=False):
    '''
    Block size directly depends on the block size of your filesystem
    to avoid performances issues
    Here I have blocks of 4096 octets (Default NTFS)
    '''
    md5 = hashlib.md5()
    with open(path,'rb') as f: 
        for chunk in iter(lambda: f.read(block_size), b''): 
             md5.update(chunk)
    if hr:
        return md5.hexdigest()

    return md5.digest()
"""
class TestCountMatrix:
    def get_count_matrix_instance(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        bucket = config.SOURCE_BUCKET
        self.key = f"{config.EXPERIMENT_ID}/r.rds"
        self.local_path = os.path.join(config.LOCAL_DIR, config.EXPERIMENT_ID)
        s3.create_bucket(
            Bucket=bucket,
            CreateBucketConfiguration={"LocationConstraint": config.AWS_REGION},
        )
        with open(os.path.join(config.LOCAL_DIR, "test", "r.rds"), "rb") as f:
            s3.upload_fileobj(Fileobj=f, Bucket=bucket, Key=self.key)

        self.count_matrix = CountMatrix()
        self.count_matrix.s3 = s3

    @mock_s3
    def test_get_objects(self):
        self.get_count_matrix_instance()
        objs = self.count_matrix.get_objects()
        assert objs == {
            "5928a56c7cbff9de78974ab50765ed20/r.rds": '"42332e7ba888546be593643616bc91fa"'
        }


    @mock_s3
    def test_download_object_does_not_exist(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = False
        Path(self.local_path).mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(self.key, "567")
        assert is_downloaded

    @mock_s3
    def test_download_object_previously_existing_etag_doesnt_match(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        Path(self.local_path).mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(self.key, "567")
        assert is_downloaded

    @mock_s3
    def test_download_object_previously_existing_etag_does_match(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        Path(self.local_path).mkdir(parents=True, exist_ok=True)
        is_downloaded = self.count_matrix.download_object(
            self.key,
            '"42332e7ba888546be593643616bc91fa"'
        )
        assert not is_downloaded

    @mock_s3
    def test_sync_no_previous_data(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True
        shutil.rmtree(self.local_path, ignore_errors=True)

        self.count_matrix.sync()
        #assert "AnnData" in type(self.count_matrix.adata).__name__
        assert not self.count_matrix.path_exists
        assert Path(self.local_path).exists()
        assert md5_for_file(
            self.adata_path
        ) == md5_for_file(
            os.path.join(config.LOCAL_DIR, "test", "r.rds")
        )

    @mock_s3
    def test_sync_previous_data_changed(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True

        shutil.rmtree(self.local_path, ignore_errors=True)
        Path(self.local_path).mkdir(parents=True, exist_ok=True)
        Path(self.adata_path).touch(exist_ok=True)

        assert md5_for_file(
                self.adata_path
        ) != md5_for_file(
                os.path.join(config.LOCAL_DIR, "test", "r.rds")
        )

        self.count_matrix.sync()
        #assert "AnnData" in type(self.count_matrix.adata).__name__
        assert self.count_matrix.path_exists
        assert Path(self.local_path).exists()
        assert md5_for_file(
            self.adata_path
        ) == md5_for_file(
            os.path.join(config.LOCAL_DIR, "test", "r.rds")
        )

    @mock_s3
    def test_sync_previous_data_not_changed(self):
        self.get_count_matrix_instance()
        self.count_matrix.path_exists = True

        Path(self.local_path).mkdir(parents=True, exist_ok=True)
        shutil.copy(
            os.path.join(config.LOCAL_DIR, "test", "r.rds"),
            self.adata_path,
        )

        self.count_matrix.sync()
        #assert "AnnData" in type(self.count_matrix.adata).__name__
        assert self.count_matrix.path_exists
        assert Path(self.local_path).exists()
        assert md5_for_file(
            self.adata_path
        ) == md5_for_file(
            os.path.join(config.LOCAL_DIR, "test", "r.rds")
        )
"""