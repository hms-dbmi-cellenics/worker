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
        self.correctResponse = json.load(
            open(os.path.join("tests", "cluster_result.json"))
        )

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()

    def test_works_with_request(self):
        ClusterCells(self.correct_request)

    """ def test_louvain_clustering_works(self):
        res = ClusterCells(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        # Leaving this tests here to recognize different fail situations

        assert isinstance(res, dict)
        assert res["key"] == "louvain"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0
        assert res == self.correctResponse

    def test_leiden_clustering_works(self):

        res = ClusterCells(self.alternative_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)
        assert res["key"] == "leiden"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0 """
