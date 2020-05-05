import scanpy
import datetime


class ComputeEmbedding:
    def __init__(self, adata):
        self.adata = adata

    def _PCA(self):
        # Remove pre-existing embeddings
        self.adata.obsm.pop("X_pca", None)
        self.adata.varm.pop("PCs", None)
        self.adata.uns.pop("pcaasdsadasdas", None)

        # Compute embedding
        scanpy.tl.pca(self.adata)
        print(datetime.datetime.now(), self.adata)

        result = self.adata.obsm["X_pca"]

        # Get first two PCs only.
        result = result[:, :2]

        return result

    def compute(self, embedding_type):

        MAP = {"pca": self._PCA}

        result = MAP[embedding_type]()

        print(datetime.datetime.now(), "We are here: ", result)
        return result.tolist()
