import json
import os

import pytest
from worker.tasks.list_genes import ListGenes


class TestListGenes:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_desc = {
            "body": {
                "name": "ListGenes",
                "selectFields": ["gene_names", "dispersions"],
                "orderBy": "dispersions",
                "orderDirection": "DESC",
                "offset": 0,
                "limit": 20,
            }
        }
        self.correct_filter = {
            "body": {
                "name": "ListGenes",
                "selectFields": ["gene_names", "dispersions"],
                "orderBy": "gene_names",
                "orderDirection": "DESC",
                "offset": 0,
                "limit": 40,
                "geneNamesFilter": "LIN",
            }
        }
        self.clean_regex = {
            "body": {
                "name": "ListGenes",
                "selectFields": ["gene_names", "dispersions"],
                "orderBy": "gene_names",
                "orderDirection": "DESC",
                "offset": 0,
                "limit": 40,
                "geneNamesFilter": "{}|()?¿*+|/.<><>LIN().?{}|()?¿*+|/.<>",
            }
        }
        self.partial_clean_regex = {
            "body": {
                "name": "ListGenes",
                "selectFields": ["gene_names", "dispersions"],
                "orderBy": "gene_names",
                "orderDirection": "DESC",
                "offset": 0,
                "limit": 40,
                "geneNamesFilter": "^${}|()?¿*+|/.<><>LIN().?{}|()?¿*+|/.<>",
            }
        }
        self.correct_response = json.load(open(os.path.join("tests", "lg_result.json")))

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ListGenes()

    def test_works_with_request(self):
        ListGenes(self.correct_desc)

    def test_construct_request(self):
        assert (
            ListGenes(self.correct_desc)._construct_request()
            == self.correct_desc["body"]
        )

    def test_construct_request_works_with_filter(self):
        assert (
            ListGenes(self.correct_filter)._construct_request()
            == self.correct_filter["body"]
        )

    def test_construct_request_cleans_regex(self):
        request = ListGenes(self.clean_regex)._construct_request()
        assert request["geneNamesFilter"] == "LIN"

    def test_construct_request_cleans_regex(self):
        request = ListGenes(self.partial_clean_regex)._construct_request()
        assert request["geneNamesFilter"] == "^$LIN"

"""
    def test_descending(self):
        res = ListGenes(self.correct_desc).compute()
        res = json.loads(res[0].result)
        assert res == self.correct_response["desc_20"]

    def test_ascending(self):
        res = ListGenes(self.correct_asc).compute()
        res = json.loads(res[0].result)
        assert res == self.correct_response["asc_20"]

    def test_by_names(self):
        res = ListGenes(self.correct_names).compute()
        res = json.loads(res[0].result)
        assert res == self.correct_response["names_desc"]

    def test_list_gene_selected_fields_appear_in_all_results(self):
        res = ListGenes(self.correct_desc).compute()
        res = res[0].result
        res = json.loads(res)

        for data in res["rows"]:
            for field in data.keys():
                assert field in self.correct_desc["body"]["selectFields"]

    def test_list_gene_has_appropriate_number_of_results(self):
        res = ListGenes(self.correct_filter).compute()
        res = res[0].result
        res = json.loads(res)
        res = res["rows"]

        assert len(res) <= self.correct_names["body"]["limit"]

    def test_filter_contains_pattern_gets_applied_to_results(self):
        res = ListGenes(self.correct_filter).compute()
        res = json.loads(res[0].result)
        for row in res["rows"]:
            assert "lin".lower() in row["gene_names"].lower()

    def test_filter_starts_with_pattern_gets_applied_to_results(self):
        res = ListGenes(self.correct_startswith).compute()
        res = json.loads(res[0].result)

        for row in res["rows"]:
            assert row["gene_names"].lower().startswith("LIN".lower())

    def test_empty_results(self):
        res = ListGenes(self.correct_empty).compute()
        res = json.loads(res[0].result)
        assert res == {"rows": [], "total": 0}
"""
