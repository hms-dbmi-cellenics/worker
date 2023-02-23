import io
import json

import boto3
import mock
import pandas as pd
import pytest
import responses
from botocore.stub import Stubber
from exceptions import RWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.cell_annotation_sctype import ScTypeAnnotate,  get_cell_sets_dict

import pdb
from worker.helpers.s3 import get_cell_sets

class TestScTypeAnnotate:
    @pytest.fixture(autouse=True)
    def get_request(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "Authorization" : "mock_authJwt",
            "body": {
                "name": "ScTypeAnnotate",
                "species": "human",
                "tissue": "Pancreas"
            }
        }

    def format_request(self):
        cell_sets = get_cell_sets(self.experiment_id)
        cell_sets_dict = get_cell_sets_dict(cell_sets)
        species = self.task_def["species"]

        return { 
            "cellSets": cell_sets_dict, 
            "species": species, 
            "tissue": tissue, 
            "apiUrl" : config.API_URL, 
            "authJwt" : "mock_authJwt"
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
        ScTypeAnnotate(self.correct_request)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ScTypeAnnotate()

    def test_generates_correct_request_keys(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = ScTypeAnnotate(self.correct_request)._format_request()
            assert isinstance(request, dict)

            # all expected keys are in the request
            expected_keys = [
                "cellSets",
                "species",
                "tissue",
                "apiUrl",
                "authJwt"
                ]
            assert all(key in request for key in expected_keys)
            stubber.assert_no_pending_responses()
 
    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/ScTypeAnnotate",
            json=payload,
            status=200,
        )

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            with pytest.raises(RWorkerException) as exception_info:
                ScTypeAnnotate(self.correct_request).compute()

            assert exception_info.value.args[0] == error_code
            assert exception_info.value.args[1] == user_message