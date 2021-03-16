import pytest
import os
import json

from tasks.cluster_cells import ClusterCells
from config import get_config

config = get_config()

# CLUSTER_ENV="test" python -m pytest --cov=. src/tests/tasks/test_cluster_cells.py


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
                "params": {},
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
                "params": {},
            },
        }
        self.correctResponse = json.load(
            open(os.path.join("tests", "cluster_result.json"))
        )

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()

    def test_louvain_clustering_works(self):
        res = ClusterCells(self.correct_request).compute()
        res = res[0].result
        print(res)
        res = json.loads(res)
        """
        I like these tests but the other one is more powerful.
        Clustering should be deterministic, let's see if this works fine, if not we can revert to old asserts. 
        assert isinstance(res, dict)
        assert res["key"] == "louvain"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0
        """
        assert res == self.correctResponse

    """
    def test_leiden_clustering_works(self):

        res = ClusterCells(self.alternative_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)
        assert res["key"] == "leiden"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0
    """
