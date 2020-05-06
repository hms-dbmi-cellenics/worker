import scanpy
import datetime
import json

from result import Result


class ComputeEmbedding:
    def __init__(self, adata):
        self.adata = adata

    def _PCA(self):
        # Remove pre-existing embeddings
        self.adata.obsm.pop("X_pca", None)
        self.adata.varm.pop("PCs", None)
        self.adata.uns.pop("pca", None)

        # Compute embedding
        scanpy.tl.pca(self.adata)
        print(datetime.datetime.now(), self.adata)

        result = self.adata.obsm["X_pca"]

        # Get first two PCs only.
        result = result[:, :2]

        return result

    def _format_result(self, result):
        # Convert numpy array to list.
        result = result.tolist()

        # JSONify list.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self, embedding_type):
        MAP = {"pca": self._PCA}

        result = MAP[embedding_type]()

        return self._format_result(result)
