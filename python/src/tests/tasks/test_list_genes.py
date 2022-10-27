import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
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

    def test_format_request(self):
        assert (
            ListGenes(self.correct_desc)._format_request()
            == self.correct_desc["body"]
        )

    def test_format_request_works_with_filter(self):
        assert (
            ListGenes(self.correct_filter)._format_request()
            == self.correct_filter["body"]
        )

    def test_format_request_cleans_regex(self):
        request = ListGenes(self.clean_regex)._format_request()
        assert request["geneNamesFilter"] == "LIN"

    def test_format_request_preserves_begin_and_end_with(self):
        request = ListGenes(self.partial_clean_regex)._format_request()
        assert request["geneNamesFilter"] == "^$LIN"
    
    @responses.activate
    def test_formats_result_appropriately(self):
        payload = {
            'data':{
                'gene_results':{
                    'gene_names': ['gene1', 'gene2', 'gene2'],
                    'dispersions': [4, 420, 1],
                },
                'full_count': 3
            }
        }

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/listGenes",
            json=payload,
            status=200,
        )
        a = ListGenes(self.correct_desc).compute().data

        pyWorkerResponse = {'total': 3, 'gene_names': ['gene1', 'gene2', 'gene2'], 'dispersions': [4, 420, 1]}
        assert (
            ListGenes(self.correct_desc).compute().data == pyWorkerResponse
        )

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/listGenes",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            ListGenes(self.correct_desc).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
