import mock
import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
from worker.helpers.mock_s3 import MockS3Class
from worker.tasks.marker_heatmap import MarkerHeatmap


class TestMarkerHeatmap:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "MarkerHeatmap",
                "cellSetKey": "set_hierarchy_1",
                "nGenes": 5,
                "type": "louvain",
                "config": {"resolution": 0.5},
            },
        }

    @pytest.fixture
    def mock_S3_get(self):
        with mock.patch("boto3.client") as m:
            mockS3 = MockS3Class()
            m.return_value = mockS3
            yield (m, mockS3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            MarkerHeatmap()

    def test_works_with_request(self):
        MarkerHeatmap(self.correct_request)

    def test_generates_correct_request_keys(self):
        MockS3Class.setResponse("hierarchichal_sets")
        request = MarkerHeatmap(self.correct_request)._format_request()
        assert isinstance(request, dict)

        # all expected keys are in the request

        expected_keys = [
            "nGenes",
            "cellSets",
        ]

        assert all(key in request for key in expected_keys)
        assert "children" in request["cellSets"].keys()
        assert request["cellSets"]["key"] == self.correct_request["body"]["cellSetKey"]

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/runMarkerHeatmap",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exc_info:
            MarkerHeatmap(self.correct_request).compute()

        assert exc_info.value.args[0] == error_code
        assert exc_info.value.args[1] == user_message
