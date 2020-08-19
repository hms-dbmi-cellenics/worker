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

        # compute embedding
        scanpy.tl.pca(self.adata)
        print(datetime.datetime.utcnow(), self.adata)
        embeddings = self.adata.obsm["X_pca"]

        # Get first two PCs only.
        raw = embeddings[:, :2]

        return raw

    def _UMAP(self):
        # Remove pre-existing embeddings
        self.adata.obsm.pop("X_umap", None)
        self.adata.uns.pop("neighbors", None)
        self.adata.obsp.pop("distances", None)
        self.adata.obsp.pop("connectivities", None)

        # Compute the neighborhood graph and create UMAP
        scanpy.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        scanpy.tl.umap(self.adata)

        print(datetime.datetime.utcnow(), self.adata)
        raw = self.adata.obsm["X_umap"]

        return raw

    def _format_result(self, raw):
        # JSONify result.
        raw_result = json.dumps(raw.tolist())

        # Return a list of formatted results.
        return [Result(raw_result), Result(raw_result)]

    def compute(self):
        MAP = {"pca": self._PCA, "umap": self._UMAP}
        embedding_type = self.task_def["type"]

        # order cells by cell ids first to guarantee order
        sorted_indices = self.adata.obs.sort_values(by=["cell_ids"]).index
        self.adata = self.adata[sorted_indices, :]

        # do the processing, get results
        raw = MAP[embedding_type]()
        return self._format_result(raw)
