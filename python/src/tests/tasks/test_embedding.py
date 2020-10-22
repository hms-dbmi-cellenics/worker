import pytest
import anndata
import os
from tasks.embedding import ComputeEmbedding
import numpy as np


class TestEmbedding:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(
            os.path.join(config.LOCAL_DIR, "test", "python.h5ad")
        )

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request_skeleton = {
            "body": {
                "name": "GetEmbedding",
                "type": "pca",
            }
        }

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            ComputeEmbedding()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            ComputeEmbedding(self.correct_request_skeleton)

    def test_works_with_request_and_adata(self):
        ComputeEmbedding(self.correct_request_skeleton, self._adata)

    def test_pca_edits_object_appropriately(self):
        try:
            old = np.array(self._adata.obsm["X_pca"][:, :2])
        except Exception:
            old = []

        res = ComputeEmbedding(self.correct_request_skeleton, self._adata)._PCA()

        assert not np.array_equal(res, old)

    def test_pca_deals_with_incomplete_previous_results(self):
        self._adata.obsm.pop("X_pca", None)
        ComputeEmbedding(self.correct_request_skeleton, self._adata)._PCA()

    def test_throws_on_invalid_embedding_type(self):
        with pytest.raises(Exception):
            ComputeEmbedding(self._adata).compute("definitelynotavalidembedding")
