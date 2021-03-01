import pytest
import anndata
import os
import json

from tasks.cluster_cells import ClusterCells
from config import get_config

config = get_config()


class TestClusterCells:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(
            os.path.join(config.LOCAL_DIR, "test", "python.h5ad")
        )

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

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            ClusterCells(self.correct_request)


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
        res = ClusterCells(alternative_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)
        assert res["key"] == "leiden"
        assert len(res["children"]) > 0
        assert len(res["children"][0]["cellIds"]) > 0
"""