import scanpy
import datetime
import json
import requests
from config import get_config
from result import Result

config = get_config()


class ComputeEmbedding:
    def __init__(self, msg, adata):
        self.adata = adata

        self.task_def = msg["body"]

    def _PCA(self):

        request = {}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getEmbeddingPCA",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        result = r.json()
        # I need to get the first two PCS of each list.
        result = [[l[0], l[1]] for l in result]
        return result

    def _UMAP(self):
        request = {}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getEmbeddingUMAP",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        result = r.json()
        return result

    def _format_result(self, raw):
        # JSONify result.
        raw_result = json.dumps(raw)

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
