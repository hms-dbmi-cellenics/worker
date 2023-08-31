import json
import pytest
import responses
from worker_status_codes import INVALID_INPUT
from exceptions import PythonWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.batch_differential_expression import BatchDifferentialExpression
from worker.helpers.s3 import get_cell_sets
from unittest.mock import patch, MagicMock
from botocore.stub import Stubber
import io
import boto3

class TestBatchDifferentialExpression:
    def get_request(
        self,
        cellSet=["cluster1"],
        compareWith="rest",
        basis=["all"],
        comparisonType=None,
    ):
        request = {
            "experimentId": config.EXPERIMENT_ID,
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "BatchDifferentialExpression",
                "cellSet": cellSet,
                "compareWith": compareWith,
                "basis": basis,
            },
        }

        if comparisonType:
            request["body"]["comparisonType"] = comparisonType

        return request

    def get_s3_stub(self, content_type):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        response = {
            "ContentLength": 10,
            "ContentType": "utf-8",
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }

        expected_params = {
            "Bucket": config.CELL_SETS_BUCKET,
            "Key": config.EXPERIMENT_ID,
        }
        print("Stubbing head_object with expected_params:")

        stubber = Stubber(s3)
        stubber.add_response("head_object", response, expected_params)

        content_bytes = json.dumps(cell_set_types[content_type], indent=2).encode(
            "utf-8"
        )
        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        response = {
            "ContentLength": len(content_bytes),
            "ContentType": "utf-8",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        return (stubber, s3)

    @responses.activate
    def test_failed_and_successful_comparison(self):
        basis = ["louvain-123", "louvain-3333"]
        stubber, s3 = self.get_s3_stub("two_sets_no_overlap")
        request_data = self.get_request(
            cellSet=["cluster1"],
            compareWith="cluster2",
            basis=basis,
            comparisonType="within",
        )

        with patch('worker.tasks.batch_differential_expression.BatchDifferentialExpression._format_request') as mock_format_request:
            valid_request = {
                "baseCells": [1, 2, 3],
                "backgroundCells": [4, 5, 6],
                "genesOnly": False,
                "comparisonType": "within",
            }
            mock_format_request.side_effect = [PythonWorkerException(INVALID_INPUT, "No data available for this comparison"), valid_request]

            with patch('requests.post') as mock_post:
                mock_post.return_value = MagicMock(status_code=200, json=lambda: {"data": {"full_count": 10, "gene_results": "Some gene results"}})

                with patch("boto3.client") as n, stubber:
                    n.return_value = s3
                    task = BatchDifferentialExpression(request_data)
                    result = task.compute()

                    assert len(result.data) == len(basis)
                    assert {"total": 0, "data": "No data available for this comparison"} in result.data
                    assert {"total": 10, "data": "Some gene results"} in result.data
