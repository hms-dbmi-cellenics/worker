import json
import os

import numpy as np
import pytest

from worker.tasks.gene_expression import GeneExpression


class TestGeneExpression:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_one_gene = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1"],
            },
        }
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "genes": ["Tpt1", "Zzz3"],
            },
        }
        self.correct_response = json.load(open(os.path.join("tests", "GE_result.json")))

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_works_with_request(self):
        GeneExpression(self.correct_request)

    def test_returns_json(self):
        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self):
        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    def test_object_returns_appropriate_number_of_genes(self):
        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == len(self.correct_request["body"]["genes"])

    def test_object_returns_one_gene(self):
        res = GeneExpression(self.correct_one_gene).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == len(self.correct_one_gene["body"]["genes"])

    def test_each_expression_data_has_correct_number_of_cells(self):
        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            assert len(v["expression"]) == 1500

    def test__expression_data_gets_displayed_appropriately(self):
        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res.values():
            expression = np.array(v["expression"], dtype=np.float)
            minimum = v["min"]
            maximum = v["max"]
            mean = v["mean"]
            stdev = v["stdev"]
            assert minimum == np.nanmin(expression)
            assert maximum == np.nanmax(expression)
            assert mean == pytest.approx(np.nanmean(expression), 0.01)
            assert stdev == pytest.approx(np.nanstd(expression), 0.01)

    # This test is commented because currently the worker doesn't handle nonexistent genes
    # A ticket has been created to fix this in expression.r
    """
    def test_task_handles_nonexistent_genes(self):

        self.correct_request["body"]["genes"] = ["PPBP", "non-existent-gene"]

        res = GeneExpression(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res) == 1
    """