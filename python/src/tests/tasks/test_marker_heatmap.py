import json
import os

import numpy as np
import pytest
from worker.tasks.marker_heatmap import MarkerHeatmap


class TestMarkerHeatmap:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "MarkerHeatmap",
                "nGenes":5,
                "type":"louvain",
                "config":{"resolution":0.5}
            },
        }
        #Correct response with flexible ngenes
        self.correct_response_genes = ['Il7r', 'S100a10', 'Crip1', 'Atp2b1', 'Il18r1', 'Cd163l1', 'Trdv4', 'Tmem176b', 'Ltb4r1', 'Tmem176a', 'Ccr7', 'Slamf6', 'Smc4', 'Gm8369', 'S100a6', 'Ccl5', 'Xcl1', 'Ly6c2', 'Klrd1', 'Nkg7', 'Gzma', 'Ncr1', 'Klra7', 'Klra9', 'Klre1', 'Cd300c2', 'Spi1', 'Alox5ap', 'Cd300a', 'Clec4a3', 'Cd79a', 'Iglc2', 'Ebf1', 'Iglc3', 'Ly6d', 'Emp2', 'Tspan7', 'Tmem100', 'Cyyr1', 'Clic5']
    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            MarkerHeatmap()

    def test_works_with_request(self):
        MarkerHeatmap(self.correct_request)

    def test_returns_json(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    """
    #Useless until we use ngenes
    def test_object_returns_appropriate_number_of_genes(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert len(res) % self.correct_request["body"]["nGenes"] == 0 
    """
    def test_object_returns_proper_genes(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        assert len(res["data"]) == len(self.correct_response_genes)

    def test_object_returns_proper_order(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        res = res["order"]
        assert res == self.correct_response_genes       

    def test_each_expression_data_has_correct_number_of_cells(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)
        res = res["data"][self.correct_response_genes[0]]["rawExpression"]
        assert len(res["expression"]) == 1500

    def test_expression_was_properly_truncated(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res["data"].values():
            truncatedExpression = np.array(v["truncatedExpression"]["expression"], dtype=np.float)
            expression = np.array(v["rawExpression"]["expression"], dtype=np.float)
            max_truncated = np.nanmax(truncatedExpression)
            lim = np.nanquantile(expression,0.95)
            i = 0.01
            while(lim==0 and i+0.95<=1):
                lim = np.nanquantile(expression,0.95+i)
                i = i+0.01

            assert max_truncated == pytest.approx(lim, 0.01)      

    def test__expression_data_gets_displayed_appropriately(self):
        res = MarkerHeatmap(self.correct_request).compute()
        res = res[0].result
        res = json.loads(res)

        for v in res["data"].values():
            expression = np.array(v["rawExpression"]["expression"], dtype=np.float)
            mean = v["rawExpression"]["mean"]
            stdev = v["rawExpression"]["stdev"]
            assert mean == pytest.approx(np.nanmean(expression), 0.01)
            assert stdev == pytest.approx(np.nanstd(expression), 0.01)
