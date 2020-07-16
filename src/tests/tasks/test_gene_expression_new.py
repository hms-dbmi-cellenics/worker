import pytest
import anndata
import os
from tasks.gene_expression_new import GeneExpressionNew
import json
from config import get_config

config = get_config()


class TestGeneExpressionNew:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(os.path.join("tests", "test.h5ad"))

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5e959f9c9f4b120771249001",
            "timeout": "2099-12-31 00:00:00",
            "body": {"name": "GeneExpressionNew", "genes": ["TGFB1", "CST3"],},
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpressionNew()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            GeneExpressionNew(self.correct_request)

    def test_works_with_request_and_adata(self):
        GeneExpressionNew(self.correct_request, self._adata)

    def test_returns_json(self):
        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self):
        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    def test_object_returns_appropriate_number_of_genes(self):
        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == len(self.correct_request["body"]["genes"])

    def test_each_expression_data_has_correct_number_of_cells(self):
        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            assert len(v["expression"]) == len(self._adata.obs)

    def test_min_and_max_expression_data_gets_displayed_appropriately(self):
        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            minimum = v["min"]
            maximum = v["max"]

            assert minimum == min(v["expression"])
            assert maximum == max(v["expression"])

    def test_task_handles_nonexistent_genes(self):

        self.correct_request["body"]["genes"] = ["TGFB1", "non-existent-gene"]

        res = GeneExpressionNew(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == 1
