import json

import mock
import pytest
import responses
from worker.config import config
from worker.tasks.differential_expression import DifferentialExpression
from worker.helpers.mock_s3 import MockS3Class

class TestDifferentialExpression:
    def get_request(
        self, cellSet="cluster1", compareWith="rest", basis="all", maxNum=None
    ):
        request = {
            "experimentId": "e52b39624588791a7889e39c617f669e1",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "DifferentialExpression",
                "cellSet": cellSet,
                "compareWith": compareWith,
                "basis": basis,
            },
        }

        if maxNum:
            request["body"]["maxNum"] = maxNum

        return request

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            json={},
            status=200,
        )

    """
    Mocks the S3 query for fetching cell sets. Returns an
    empty cell set and yields the patched up object.
    """

    @pytest.fixture
    def mock_S3_get(self):
        with mock.patch("boto3.client") as m:
            mockS3 = MockS3Class()
            m.return_value = mockS3
            yield (m, mockS3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    @responses.activate
    def test_generates_correct_request_keys(self, mock_S3_get):
        MockS3Class.setResponse("one_set")
        request = DifferentialExpression(self.get_request())._format_request()
        print(request)
        assert isinstance(request, dict)

        # all expected keys are in the request
        expected_keys = [
            "baseCells",
            "backgroundCells",
        ]
        
        assert all(key in request for key in expected_keys)
