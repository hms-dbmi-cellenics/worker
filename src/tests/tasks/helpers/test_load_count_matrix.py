from helpers.load_count_matrix import _load_file
from botocore.stub import Stubber, ANY
from moto import mock_s3
import boto3
import mock
from config import get_config

from helpers.dynamo import get_item_from_dynamo

config = get_config()


class TestLoadCountMatrix:
    @mock_s3
    def test_returns_correct_adata_object_when_path_and_key_exist_in_non_development(
        self,
    ):
        s3 = boto3.client("s3")
        bucket = "my_custom_bucket_path"
        key = "very/long/and/convoluted/path"
        s3.create_bucket(Bucket=bucket)

        with open("tests/test.h5ad", "rb") as f, mock.patch("config.get_config") as m:
            mock_config = config
            mock_config.ENVIRONMENT = "staging"
            m.return_value = mock_config

            s3.upload_fileobj(f, bucket, key)
            a = _load_file(f"{bucket}/{key}")

            assert "AnnData" in type(a).__name__
