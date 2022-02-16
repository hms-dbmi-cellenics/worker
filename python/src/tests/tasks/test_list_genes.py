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

    def test_construct_request_preserves_begin_and_end_with(self):
        request = ListGenes(self.partial_clean_regex)._construct_request()
        assert request["geneNamesFilter"] == "^$LIN"
