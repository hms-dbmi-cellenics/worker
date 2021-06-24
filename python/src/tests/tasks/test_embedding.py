import pytest
import os
from worker.tasks.embedding import GetEmbedding
import json


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
                "config": {"minimumDistance": 0.5, "distanceMetric": "euclidean"},
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
        """
        The test file has been created with the multisample dataset, expId: e52b39624588791a7889e39c617f669e
        """
        self.correctResponse = json.load(open(os.path.join("tests", "emb_result.json")))

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GetEmbedding()

    def test_works_with_request(self):
        GetEmbedding(self.correct_request_pca)

    def test_pca_edits_object_appropriately(self):

        old = str(self.correctResponse["pca"])
        res = GetEmbedding(self.correct_request_pca).compute()[0].result

        assert res == old

    def test_umap_edits_object_appropriately(self):

        old = str(self.correctResponse["umap"])
        res = GetEmbedding(self.correct_request_umap).compute()[0].result

        assert res == old

    def test_umap_different_params(self):

        old = str(self.correctResponse["umap_cosine"])
        res = GetEmbedding(self.correct_request_umap_cosine).compute()[0].result

        assert res == old

    def test_tsne_edits_object_appropriately(self):

        old = str(self.correctResponse["tsne"])
        res = GetEmbedding(self.correct_request_tsne).compute()[0].result

        assert res == old

    def test_throws_on_invalid_task_def(self):
        with pytest.raises(Exception):
            GetEmbedding(self).compute("definitelynotavalidembedding")
