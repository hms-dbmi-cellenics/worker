import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.cluster_cells import ClusterCells


class TestClusterCells:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "Authorization" : "mock_authJwt",
            "body": {
                "name": "ClusterCells",
                "cellSetName": "Louvain clusters",
                "type": "louvain",
                "cellSetKey": "louvain",
                "config": {"resolution": 0.5},
            },
        }
        self.parsed_request = {
            "type": self.correct_request["body"]["type"],
            "config": self.correct_request["body"]["config"],
            "apiUrl": config.API_URL,
            "authJwt": "mock_authJwt"
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ClusterCells()

    def test_works_with_request(self):
        ClusterCells(self.correct_request)

    def test_format_request(self):
        assert (
            ClusterCells(self.correct_request)._format_request() == self.parsed_request
        )

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/getClusters",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            ClusterCells(self.correct_request).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
