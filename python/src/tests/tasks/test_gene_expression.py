import json
import io

import random
import numpy as np
import boto3
import pytest
import responses
import mock

from botocore.stub import Stubber
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.gene_expression import GeneExpression

from tests.data.cell_sets_from_s3 import cell_sets_from_s3

from tests.utils import get_cell_ids 

class TestGeneExpression:
    def get_s3_stub(self, cell_sets):
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
        content_bytes = json.dumps(cell_sets, indent=2).encode(
            "utf-8"
        )

        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        response = {
            "ContentLength": len(cell_sets),
            "ContentType": "utf-8",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        return (stubber, s3)

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_one_gene = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
            },
        }
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1", "Zzz3"],
            },
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_works_with_request(self):
        GeneExpression(self.correct_request)

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/runExpression",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            GeneExpression(self.correct_request).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
