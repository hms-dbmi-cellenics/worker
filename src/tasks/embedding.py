import scanpy
import datetime
import json

from result import Result


class ComputeEmbedding:
    def __init__(self, msg, adata):
        self.adata = adata

        self.task_def = msg["body"]

    def _PCA(self):
        # Remove pre-existing embeddings
        self.adata.obsm.pop("X_pca", None)
        self.adata.varm.pop("PCs", None)
        self.adata.uns.pop("pca", None)

        # Compute embedding
        scanpy.tl.pca(self.adata)
        print(datetime.datetime.now(), self.adata)
        embeddings = self.adata.obsm["X_pca"]

        # Get first two PCs only.
        embeddings = embeddings[:, :2]

        result = {}

        for index, data in zip(self.adata.obs.index, embeddings):
            result[index] = data.tolist()

        return result

    def _format_result(self, result):

        print("we are dumping", result)

        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        MAP = {"pca": self._PCA}

        embedding_type = self.task_def["type"]
        result = MAP[embedding_type]()
        return self._format_result(result)
