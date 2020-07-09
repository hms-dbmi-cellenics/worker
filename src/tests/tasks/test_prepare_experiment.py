import pytest
from botocore.stub import Stubber, ANY

import boto3
import mock
from config import get_config
import shutil
import os

from tasks.prepare_experiment import PrepareExperiment


config = get_config()


class TestPrepareExperiment:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request_skeleton = {
            "experimentId": "12345",
            "body": {
                "name": "PrepareExperiment",
                "sourceBucket": "biomage-source-originals",
                "sourceMatrixPath": "fake/",
            },
        }

    def test_prepare_experiment(self):
        path = self.correct_request_skeleton["body"]["sourceMatrixPath"]
        shutil.copytree("tests/count_mtrx", path)
        adata = None
        instance = PrepareExperiment(self.correct_request_skeleton, adata)

        s3 = boto3.client("s3")
        stubber = Stubber(s3)
        stubber.add_response(
            "list_objects",
            {
                "Contents": [
                    {"Key": path, "Size": 0},
                    {"Key": "/".join([path, "barcodes.tsv"]), "Size": 12533760},
                    {"Key": "/".join([path, "genes.tsv"]), "Size": 816952},
                    {"Key": "/".join([path, "matrix.mtx"]), "Size": 42362695},
                ]
            },
            {
                "Bucket": self.correct_request_skeleton["body"]["sourceBucket"],
                "Prefix": path,
            },
        )

        s3.download_fileobj = mock.MagicMock()
        s3.upload_fileobj = mock.MagicMock()
        stubber.activate()

        dynamo = boto3.resource("dynamodb")
        dynamo_stubber = Stubber(dynamo.meta.client)
        dynamo_stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.get_dynamo_table(),
                "Key": {"experimentId": self.correct_request_skeleton["experimentId"]},
                "ProjectionExpression": "matrixPath",
            },
        )
        dynamo_stubber.activate()

        with mock.patch("boto3.client") as m, mock.patch("boto3.resource") as m2:
            m.return_value = s3
            m2.return_value = dynamo
            r = instance.compute()
            assert r[0].result == '{"url": "https://my-experiment-url/12345"}'

        assert not os.path.exists(path)

    def test_prepare_experiment_fails_dynamo_lookup(self):
        path = self.correct_request_skeleton["body"]["sourceMatrixPath"]
        shutil.copytree("tests/count_mtrx", path)
        adata = None
        instance = PrepareExperiment(self.correct_request_skeleton, adata)

        s3 = boto3.client("s3")
        stubber = Stubber(s3)
        stubber.add_response(
            "list_objects",
            {
                "Contents": [
                    {"Key": path, "Size": 0},
                    {"Key": "/".join([path, "barcodes.tsv"]), "Size": 12533760},
                    {"Key": "/".join([path, "genes.tsv"]), "Size": 816952},
                    {"Key": "/".join([path, "matrix.mtx"]), "Size": 42362695},
                ]
            },
            {
                "Bucket": self.correct_request_skeleton["body"]["sourceBucket"],
                "Prefix": path,
            },
        )

        s3.download_fileobj = mock.MagicMock()
        s3.upload_fileobj = mock.MagicMock()
        stubber.activate()

        dynamo = boto3.resource("dynamodb")
        dynamo_stubber = Stubber(dynamo.meta.client)
        dynamo_stubber.add_response(
            "get_item",
            {},
            {
                "TableName": config.get_dynamo_table(),
                "Key": {"experimentId": self.correct_request_skeleton["experimentId"]},
                "ProjectionExpression": "matrixPath",
            },
        )
        dynamo_stubber.activate()

        with mock.patch("boto3.client") as m, mock.patch("boto3.resource") as m2:
            m.return_value = s3
            m2.return_value = dynamo
            with pytest.raises(Exception):
                instance.compute()

        assert not os.path.exists(path)
