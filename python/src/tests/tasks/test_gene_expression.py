import pytest
import anndata
import os
import statistics
from tasks.gene_expression import GeneExpression
import json
from config import get_config

config = get_config()


class TestGeneExpression:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(
            os.path.join(config.LOCAL_DIR, "test", "python.h5ad")
        )

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5e959f9c9f4b120771249001",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["PPBP", "CST3"],
            },
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            GeneExpression(self.correct_request)

    def test_works_with_request_and_adata(self):
        GeneExpression(self.correct_request, self._adata)

    def test_returns_json(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    def test_object_returns_appropriate_number_of_genes(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == len(self.correct_request["body"]["genes"])

    def test_each_expression_data_has_correct_number_of_cells(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            assert len(v["expression"]) == len(self._adata.obs)

    def test__expression_data_gets_displayed_appropriately(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            minimum = v["min"]
            maximum = v["max"]
            mean = v["mean"]
            stdev = v["stdev"]

            assert minimum == min(v["expression"])
            assert maximum == max(v["expression"])
            assert mean == statistics.mean(v["expression"])
            assert stdev == statistics.stdev(v["expression"])

    def test_task_handles_nonexistent_genes(self):

        self.correct_request["body"]["genes"] = ["PPBP", "non-existent-gene"]

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == 1
