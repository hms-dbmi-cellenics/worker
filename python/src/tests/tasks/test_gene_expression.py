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
    def test_should_not_parse_json_error_in_200_response(self):
        # GeneExpression intentionally doesn't parse the full JSON response
        # to avoid parsing very large response bodies
        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/runExpression",
            json=payload,
            status=200,
        )

        # Should not raise an exception since we don't parse the JSON
        result = GeneExpression(self.correct_request).compute()
        assert result is not None

    @responses.activate
    def test_should_throw_exception_on_error_status(self):
        # Verify that HTTP error status codes raise an exception
        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/runExpression",
            json={"error": "Server error"},
            status=500,
        )

        with pytest.raises(Exception):  # raise_for_status raises HTTPError
            GeneExpression(self.correct_request).compute()

    def test_format_request_sets_cell_ids_to_none_without_downsample_settings(
        self,
    ):
        task = GeneExpression(self.correct_request)
        request = task._format_request()
        assert request["cellIds"] is None

    def test_format_request_uses_provided_cell_ids(self):
        msg = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
                "downsampleSettings": {
                    "cellIds": [1, 2, 3],
                    "selectedCellSet": "louvain",
                    "groupedTracks": ["sample"],
                },
            },
        }
        task = GeneExpression(msg)
        request = task._format_request()
        assert request["cellIds"] == [1, 2, 3]

    def test_format_request_fetches_bucketed_cells_when_no_cell_ids(self):
        msg = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
                "downsampleSettings": {
                    "selectedCellSet": "louvain",
                    "groupedTracks": ["sample"],
                },
            },
        }
        bucketed_ids = [10, 20, 30]
        with mock.patch(
            "worker.tasks.gene_expression.get_cell_sets",
            return_value=cell_sets_from_s3["cellSets"],
        ), mock.patch(
            "worker.tasks.gene_expression.get_bucketed_heatmap_cells",
            return_value=bucketed_ids,
        ) as mock_bucketed:
            task = GeneExpression(msg)
            request = task._format_request()
            mock_bucketed.assert_called_once_with(
                "louvain",
                ["sample"],
                cell_sets_from_s3["cellSets"],
            )
            assert request["cellIds"] == bucketed_ids
