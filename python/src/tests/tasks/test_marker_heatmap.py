import io
import json

import boto3
import mock
import pytest
import responses
from botocore.stub import Stubber
from exceptions import RWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.marker_heatmap import MarkerHeatmap


class TestMarkerHeatmap:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": config.EXPERIMENT_ID,
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "MarkerHeatmap",
                "cellSetKey": "set_hierarchy_1",
                "nGenes": 5,
                "type": "louvain",
                "config": {"resolution": 0.5},
            },
        }

    """
    Returns a stubber and a stubbed s3 client that will get executed
    in the code instead of the real s3 clients and return the desired
    cell sets content, depending on content_type
    """

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
            "ContentLength": len(cell_set_types[content_type]),
            "ContentType": "utf-8",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        return (stubber, s3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            MarkerHeatmap()

    def test_works_with_request(self):
        MarkerHeatmap(self.correct_request)

    def test_generates_correct_request_keys(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            bla = MarkerHeatmap(self.correct_request)

            request = bla._format_request()
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
