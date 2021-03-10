import pytest
import os
import statistics
from tasks.doublet_score import GetDoubletScore
import json
from config import get_config
import responses

config = get_config()


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

    @pytest.fixture(autouse=True)
    def set_responses(self):
        with open(os.path.join("tests", "DoubletScore_result.json")) as f:
            data = json.load(f)
            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/getDoubletScore",
                json=data,
                status=200,
            )


    def test_works_with_request(self):
        GetDoubletScore()

    @responses.activate
    def test_returns_json(self):
        res = GetDoubletScore().compute()
        res = res[0].result
        json.loads(res)

    @responses.activate
    def test_returns_a_json_object(self):
        res = GetDoubletScore().compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, list)