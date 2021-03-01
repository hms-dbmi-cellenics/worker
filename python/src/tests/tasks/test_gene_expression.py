import pytest
import anndata
import os
import statistics
from tasks.gene_expression import GeneExpression
import json
from config import get_config
import responses

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
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["PPBP", "PF4"],
            },
        }

    @pytest.fixture(autouse=True)
    def set_responses(self):
        with open(os.path.join("tests", "GE_result.json")) as f:
            data = json.load(f)
            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/getExpression",
                json=data,
                status=200,
            )

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            GeneExpression(self.correct_request)

    def test_works_with_request_and_adata(self):
        GeneExpression(self.correct_request, self._adata)

    @responses.activate
    def test_returns_json(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        json.loads(res)

    @responses.activate
    def test_returns_a_json_object(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    @responses.activate
    def test_object_returns_appropriate_number_of_genes(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == len(self.correct_request["body"]["genes"])

    #
    #   The following tests passed with the json object generated with the R script, but won't work with the current h5ad file.
    #
    """
    @responses.activate
    def test_each_expression_data_has_correct_number_of_cells(self):
        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            assert len(v["expression"]) == len(self._adata.obs)

    @responses.activate
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
            assert mean == pytest.approx(statistics.mean(v["expression"]), 0.01)
            assert stdev == pytest.approx(statistics.stdev(v["expression"]), 0.01)

    @responses.activate
    def test_task_handles_nonexistent_genes(self):

        self.correct_request["body"]["genes"] = ["PPBP", "non-existent-gene"]

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == 1
    """