import pytest
import anndata
import os
import json
from tasks.list_genes import ListGenes
from config import get_config

config = get_config()


class TestListGenes:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(
            os.path.join(config.LOCAL_DIR, "test", "python.h5ad")
        )

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request_skeleton = {
            "body": {
                "name": "ListGenes",
                "selectFields": ["gene_names", "dispersions"],
                "orderBy": "dispersions",
                "orderDirection": "desc",
                "offset": 0,
                "limit": 20,
            }
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ListGenes()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            ListGenes(self.correct_request_skeleton)

    def test_works_with_request_and_adata(self):
        ListGenes(self.correct_request_skeleton, self._adata)

    def test_returns_json(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert isinstance(res, dict)

    def test_result_object_has_total_which_is_int(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert isinstance(res["total"], int)

    def test_result_object_has_total_results_which_is_list(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert isinstance(res["rows"], list)

    def test_list_gene_selected_fields_appear_in_all_results(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for data in res["rows"]:
            for field in data.keys():
                assert field in self.correct_request_skeleton["body"]["selectFields"]

    def test_list_gene_has_appropriate_number_of_results(self):
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        res = res["rows"]

        assert len(res) <= self.correct_request_skeleton["body"]["limit"]

    def test_filter_contains_pattern_gets_applied_to_results(self):
        task_def = self.correct_request_skeleton["body"]
        task_def["geneNamesFilter"] = "%LIN%"
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = json.loads(res[0].result)

        for row in res["rows"]:
            assert "LIN" in row["gene_names"]

    def test_filter_starts_with_pattern_gets_applied_to_results(self):
        task_def = self.correct_request_skeleton["body"]
        task_def["geneNamesFilter"] = "LIN%"
        res = ListGenes(self.correct_request_skeleton, self._adata).compute()
        res = json.loads(res[0].result)

        for row in res["rows"]:
            assert row["gene_names"].startswith("LIN")
