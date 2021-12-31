import json
import os

import mock
import numpy as np
import pytest
import responses
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
        # Correct response with flexible ngenes
        self.correct_response_genes = [
            "Il7r",
            "S100a10",
            "Crip1",
            "Atp2b1",
            "Il18r1",
            "Cd163l1",
            "Trdv4",
            "Tmem176b",
            "Ltb4r1",
            "Tmem176a",
            "Ccr7",
            "Slamf6",
            "Smc4",
            "Gm8369",
            "S100a6",
            "Ccl5",
            "Xcl1",
            "Ly6c2",
            "Klrd1",
            "Nkg7",
            "Gzma",
            "Ncr1",
            "Klra7",
            "Klra9",
            "Klre1",
            "Cd300c2",
            "Spi1",
            "Alox5ap",
            "Cd300a",
            "Clec4a3",
            "Cd79a",
            "Iglc2",
            "Ebf1",
            "Iglc3",
            "Ly6d",
            "Emp2",
            "Tspan7",
            "Tmem100",
            "Cyyr1",
            "Clic5",
        ]

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

    @responses.activate
    def test_generates_correct_request_keys(self, mock_S3_get):
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
