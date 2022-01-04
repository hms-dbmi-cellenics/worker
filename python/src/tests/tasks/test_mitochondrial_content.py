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
            }
        }

    def test_works_with_request(self):
        GetMitochondrialContent(self.correct_request)

    def test_generates_correct_request_keys(self):
        request = GetDoubletScore(self.correct_request)._format_request()
        assert isinstance(request, dict)

        # all expected keys are in the request
        expected_keys = []
        assert all(key in request for key in expected_keys)

