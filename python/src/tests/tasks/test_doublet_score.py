import json
import os

import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.doublet_score import GetDoubletScore


class TestGetDoubletScore:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5e959f9c9f4b120771249001",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "getDoubletScore",
            },
        }

    def test_works_with_request(self):
        GetDoubletScore(self.correct_request)

    def test_generates_correct_request_keys(self):
        request = GetDoubletScore(self.correct_request)._format_request()
        assert isinstance(request, dict)

        # all expected keys are in the request
        expected_keys = []
        assert all(key in request for key in expected_keys)

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/getDoubletScore",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            GetDoubletScore(self.correct_request).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
