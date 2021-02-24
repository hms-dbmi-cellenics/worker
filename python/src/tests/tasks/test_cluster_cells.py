import pytest
import os
import json

from tasks.cluster_cells import ClusterCells
from config import get_config

config = get_config()


class TestClusterCells:

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5928a56c7cbff9de78974ab50765ed20",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "ClusterCells",
                "cellSetName": "Louvain clusters",
                "type": "louvain",
                "cellSetKey": "louvain",
                "params": {},
            },
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()


"""
    def test_louvain_clustering_works(self):
        res = ClusterCells(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)
        assert res["key"] == "louvain"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0

    def test_leiden_clustering_works(self):
        alternative_request = {
            "experimentId": "5928a56c7cbff9de78974ab50765ed20",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "ClusterCells",
                "cellSetName": "Leiden clusters",
                "type": "leiden",
                "cellSetKey": "leiden",
                "params": {},
            },
        }
        res = ClusterCells(alternative_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)
        assert res["key"] == "leiden"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0
"""