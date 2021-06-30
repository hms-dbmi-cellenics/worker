import json
import os

import pytest
import responses

from worker.config import config
from worker.tasks.mitochondrial_content import GetMitochondrialContent


class TestGetMitochondrialContent:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5e959f9c9f4b120771249001",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "getMitochondrialContent",
            },
        }

    @pytest.fixture(autouse=True)
    def set_responses(self):
        with open(os.path.join("tests", "Mitochondrial_result.json")) as f:
            data = json.load(f)
            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
                json=data,
                status=200,
            )

    def test_works_with_request(self):
        GetMitochondrialContent(self.correct_request)

    @responses.activate
    def test_returns_json(self):
        res = GetMitochondrialContent(self.correct_request).compute()
        res = res[0].result
        json.loads(res)

    @responses.activate
    def test_returns_a_json_object(self):
        res = GetMitochondrialContent(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, list)
