from moto import mock_s3
import boto3
from helpers.load_count_matrix import _load_file
from config import get_config

config = get_config()


class TestLoadCountMatrix:
    @mock_s3
    def test_returns_correct_adata_object_when_path_and_key_exist(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        bucket = "my-custom-bucket-path"
        key = "very/long/and/convoluted/path"
        s3.create_bucket(
            Bucket=bucket,
            CreateBucketConfiguration={"LocationConstraint": config.AWS_REGION},
        )

        with open("tests/test.h5ad", "rb") as f:
            s3.upload_fileobj(Fileobj=f, Bucket=bucket, Key=key)
            a = _load_file(f"{bucket}/{key}")

        assert "AnnData" in type(a).__name__
