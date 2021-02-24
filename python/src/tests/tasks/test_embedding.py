import pytest
import anndata
import os
import numpy as np
from tasks.embedding import ComputeEmbedding
from config import get_config
import json
import responses

config = get_config()


class TestEmbedding:

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request_skeleton = {
            "body": {"name": "GetEmbedding", "type": "pca", "config": ""}
        }

    @pytest.fixture(autouse=True)
    def set_responses(self):
        with open(os.path.join("tests", "emb_result.json")) as f:
            data = json.load(f)
            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/getEmbedding",
                json=data,
                status=200,
            )

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ComputeEmbedding()

    def test_works_with_request(self):
        ComputeEmbedding(self.correct_request_skeleton)

    #
    # These two tests are not useful in R. Will replace with additional testing when we pick the correct r testing framework.
    #
    """
    @responses.activate
    def test_pca_edits_object_appropriately(self):
        try:
            old = np.array(self._adata.obsm["X_pca"][:, :2])
        except Exception:
            old = []

        res = ComputeEmbedding(self.correct_request_skeleton, self._adata).compute()

        assert True

    @responses.activate
    def test_pca_deals_with_incomplete_previous_results(self):
        self._adata.obsm.pop("X_pca", None)
        ComputeEmbedding(self.correct_request_skeleton, self._adata).compute()

    def test_throws_on_invalid_embedding_type(self):
        with pytest.raises(Exception):
            ComputeEmbedding(self._adata).compute("definitelynotavalidembedding")
"""