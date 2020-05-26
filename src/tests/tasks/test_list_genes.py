import pytest
import anndata
import os
from tasks.list_genes import ListGenes
from result import Result
import numpy as np
import json


class TestEmbedding:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(os.path.join("tests", "test.h5ad"))

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_task_def = {
            "name": "ListGenes",
            "selectFields": ["highly_variable", "dispersions"],
            "orderBy": "dispersions",
            "orderDirection": "desc",
            "offset": 0,
            "limit": 20,
        }

    def test_list_genes_throws_on_missing_anndata(self):
        with pytest.raises(TypeError):
            ListGenes()

    def test_list_genes_works_with_test_data(self):
        ListGenes(self._adata)

    def test_list_genes_returns_json(self):
        res = ListGenes(self._adata).compute(self.correct_task_def)
        res = res[0].result
        json.loads(res)

    def test_list_genes_returns_a_json_list(self):
        res = ListGenes(self._adata).compute(self.correct_task_def)
        res = res[0].result
        res = json.loads(res)

        assert isinstance(res, list)

    def test_list_gene_selected_fiels_appear_in_all_results(self):
        res = ListGenes(self._adata).compute(self.correct_task_def)
        res = res[0].result
        res = json.loads(res)

        for data in res:
            for field in data.keys():
                assert field in self.correct_task_def["selectFields"]

    def test_list_gene_has_appropriate_number_of_results(self):
        res = ListGenes(self._adata).compute(self.correct_task_def)
        res = res[0].result
        res = json.loads(res)

        assert len(res) <= self.correct_task_def["limit"]
