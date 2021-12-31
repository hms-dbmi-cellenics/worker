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
