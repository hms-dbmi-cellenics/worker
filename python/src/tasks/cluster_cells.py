import scanpy
import json

from result import Result
from helpers.color_pool import COLOR_POOL


class ClusterCells:
    def __init__(self, msg, adata):
        self.adata = adata

        self.experiment_id = msg["experimentId"]
        self.task_def = msg["body"]
        self.colors = COLOR_POOL.copy()

    def _louvain(self, params):
        scanpy.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        scanpy.tl.louvain(self.adata, key_added="cluster", **params)

        return self.adata.obs[["cluster", "cell_ids"]]

    def _leiden(self, params):
        scanpy.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        scanpy.tl.leiden(self.adata, key_added="cluster", **params)

        return self.adata.obs[["cluster", "cell_ids"]]

    def _format_result(self, raw, cell_set_key, cell_set_name):
        # construct new cell set group
        cell_set = {
            "key": cell_set_key,
            "name": cell_set_name,
            "rootNode": True,
            "children": [],
        }

        for cluster in raw["cluster"].cat.categories:
            view = raw[raw.cluster == cluster]["cell_ids"]
            cell_set["children"].append(
                {
                    "key": f"{cell_set_key}-{cluster}",
                    "name": f"Cluster {cluster}",
                    "color": self.colors.pop(0),
                    "cellIds": [int(id) for id in view.tolist()],
                }
            )

        return [Result(json.dumps(cell_set), cacheable=False)]

    def compute(self):
        cluster_type = self.task_def["type"]
        cell_set_key = self.task_def["cellSetKey"]
        cell_set_name = self.task_def["cellSetName"]

        MAP = {"leiden": self._leiden, "louvain": self._louvain}
        params = self.task_def.get("params", {})

        # do the processing, get results
        raw = MAP[cluster_type](params)
        return self._format_result(raw, cell_set_key, cell_set_name)
