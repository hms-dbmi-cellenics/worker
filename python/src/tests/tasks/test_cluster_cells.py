import json
import os

import pytest
from worker.tasks.cluster_cells import ClusterCells


class TestClusterCells:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "ClusterCells",
                "cellSetName": "Louvain clusters",
                "type": "louvain",
                "cellSetKey": "louvain",
                "config": {"resolution": 0.5},
            },
        }
        self.alternative_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "ClusterCells",
                "cellSetName": "Leiden clusters",
                "type": "leiden",
                "cellSetKey": "leiden",
                "config": {"resolution": 0.5},
            },
        }
        self.parsed_request = {
            "type": self.correct_request["body"]["type"],
            "config": self.correct_request["body"]["config"],
        }

        self.correctResponse = json.load(
            open(os.path.join("tests", "cluster_result.json"))
        )

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()

    def test_works_with_request(self):
        ClusterCells(self.correct_request)

    def test_format_request(self):
        assert (
            ClusterCells(self.correct_request)._format_request()
            == self.parsed_request
        )
