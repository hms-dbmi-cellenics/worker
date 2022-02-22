import io
import json
import os

import boto3
import mock
import pytest
from botocore.stub import Stubber

from worker.config import config
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
        with open(os.path.join("tests/data", "cell_set_types.json")) as f:
            all_cell_set_types = json.load(f)
        content = all_cell_set_types[content_type]
        content_bytes = json.dumps(content, indent=2).encode("utf-8")
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

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            MarkerHeatmap()

    def test_works_with_request(self):
        MarkerHeatmap(self.correct_request)

    def test_generates_correct_request_keys(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = MarkerHeatmap(self.correct_request)._format_request()
            assert isinstance(request, dict)

            # all expected keys are in the request

            expected_keys = [
                "nGenes",
                "cellSets",
            ]

            assert all(key in request for key in expected_keys)
            assert "children" in request["cellSets"].keys()
            assert (
                request["cellSets"]["key"] == self.correct_request["body"]["cellSetKey"]
            )
            stubber.assert_no_pending_responses()
