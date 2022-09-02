import io
import json

import boto3
import mock
import pytest
import responses
import pandas as pd
from botocore.stub import Stubber

from exceptions import RWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.normalized_matrix import GetNormalizedExpression
from worker.tasks.dotplot import DotPlot


class TestGetNormalizedExpression:
    @pytest.fixture(autouse=True)
    def get_request(self):
        self.correct_request = {
            "body": {
                "name": "GetNormalizedExpression",
                "filterBy":{
                    "group":["set_hierarchy_2"],
                    "key":["cluster4"]
                }
            }
        }


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

        # Get object
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

    def test_works_with_request(self):
        GetNormalizedExpression(self.correct_request)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GetNormalizedExpression()

    def test_generates_correct_request_keys(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = GetNormalizedExpression(self.correct_request)._format_request()
            assert isinstance(request, dict)

            # all expected keys are in the request
            expected_keys = [
                "filterBy",
                "applyFilter",
                ]
            assert all(key in request for key in expected_keys)
            stubber.assert_no_pending_responses()

    def test_generates_correct_number_of_cellids(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = GetNormalizedExpression(self.correct_request)._format_request()
             # Get object
            cell_set = cell_set_types["hierarchichal_sets"]

            assert len(request["filterBy"]) == len(cell_set["cellSets"][1]["children"][1]["cellIds"])
            stubber.assert_no_pending_responses()

    @responses.activate
    def test_generates_correct_result_type(self):
        payload = {'data': {'CACGGGTTCTGTTGGA-1': [0,0], 'CACGGGTTCTGTTGGT-1': [0,0],
        'CACGGGTTCTGTTGTA-1': [0,0], 'CACGGGTTCTGTTTGA-1': [0,0],
        'CACGGGTTCTGTTGGC-1': [0,0], 'CACGGGTTCTGTTGCC-1': [0,0],
        '_row':['ENSG0', 'ENSG1']}}

        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3

            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/GetNormalizedExpression",
                json=payload,
                status=200,
            )

            result = GetNormalizedExpression(self.correct_request).compute()
            assert isinstance(result.data, pd.DataFrame)
            stubber.assert_no_pending_responses()

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/GetNormalizedExpression",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            GetNormalizedExpression(self.correct_request).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
