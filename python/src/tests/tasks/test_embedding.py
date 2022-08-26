import pytest
import responses
from exceptions import RWorkerException
from worker.config import config
from worker.tasks.embedding import GetEmbedding


class TestEmbedding:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request_pca = {
            "body": {"name": "GetEmbedding", "type": "pca", "config": ""}
        }
        self.correct_request_umap = {
            "body": {
                "name": "GetEmbedding",
                "type": "umap",
                "config": {
                    "minimumDistance": 0.5,
                    "distanceMetric": "euclidean",
                },
            }
        }
        self.correct_request_umap_cosine = {
            "body": {
                "name": "GetEmbedding",
                "type": "umap",
                "config": {"minimumDistance": 0.1, "distanceMetric": "cosine"},
            }
        }
        self.correct_request_tsne = {
            "body": {
                "name": "GetEmbedding",
                "type": "tsne",
                "config": {"perplexity": 30, "learningRate": 200},
            }
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GetEmbedding()

    def test_works_with_request(self):
        GetEmbedding(self.correct_request_pca)

    def test_throws_on_invalid_task_def(self):
        with pytest.raises(Exception):
            GetEmbedding(self).compute("definitelynotavalidembedding")

    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

        error_code = "MOCK_R_WORKER_ERROR"
        user_message = "Some worker error"

        payload = {"error": {"error_code": error_code, "user_message": user_message}}

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/getEmbedding",
            json=payload,
            status=200,
        )

        with pytest.raises(RWorkerException) as exception_info:
            GetEmbedding(self.correct_request_umap).compute()

        assert exception_info.value.args[0] == error_code
        assert exception_info.value.args[1] == user_message
