import json
import os

import numpy as np
import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.gene_expression import GeneExpression


class TestGeneExpression:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_one_gene = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
                "downsampled": False,
            },
        }
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1", "Zzz3"],
                "downsampled": False,
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
