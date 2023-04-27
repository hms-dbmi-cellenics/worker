import json
import pytest
import responses
from exceptions import ErrorCodes, PythonWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.batch_differential_expression import BatchDifferentialExpression
from worker.helpers.s3 import get_cell_sets
from unittest.mock import patch
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
    
    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            BatchDifferentialExpression()
    
    @responses.activate
    def test_should_handle_r_worker_error(self):
        basis = ["louvain-123", "louvain-3333"]
        stubber, s3 = self.get_s3_stub("two_sets_no_overlap")
        error_code = "R_WORKER_ERROR"
        user_message = "User message"
        request_data = self.get_request(
            cellSet=["cluster1"],
            compareWith="cluster2",
            basis=basis,
            comparisonType="within",
        )

        with patch('worker.helpers.get_diff_expr_cellsets') as mock_func:
            # Create a list of exceptions with the same length as the number of basis values
            exceptions = [PythonWorkerException(ErrorCodes.INVALID_INPUT, "No cell id fullfills the 1st cell set.") for _ in basis]
            mock_func.side_effect = exceptions

            with patch("boto3.client") as n, stubber:
                n.return_value = s3
                task = BatchDifferentialExpression(request_data)
                result = task.compute()
                print('RESULT IS ', result)
                assert {"total": 0, "data": "No cell id fulfills the 1st cell set."} in result.data



    # def test_should_call_r_worker_endpoint(self):
    #     cell_sets = ["louvain-0","louvain-1","louvain-2","louvain-3","louvain-4"]
    #     request_data = self.get_request(
    #         cellSet=cell_sets,
    #         compareWith="background",
    #         basis=["all"],
    #         comparisonType="within",
    #     )

    #     with patch("worker.helpers.s3.get_cell_sets", return_value=cell_set_types["three_sets"]), \
    #         responses.RequestsMock() as rsps:
    #             # Add mocked requests for each iteration in the loop
    #             for _ in range(len(cell_sets)):  # Assuming 4 requests will be made based on the 'basis' value
    #                 rsps.add(
    #                     responses.POST,
    #                     f"{config.R_WORKER_URL}/v0/DifferentialExpression",
    #                     json={"data": {"full_count": 10, "gene_results": ['ok']}},
    #                     status=200,
    #                 )

    #             task = BatchDifferentialExpression(request_data)
    #             result = task.compute()
    #             assert {"total": 10, "data": ['ok']} in result.data