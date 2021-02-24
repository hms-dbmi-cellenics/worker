import pandas as pd
import json
import requests
from result import Result
from helpers.color_pool import COLOR_POOL
from config import get_config

config = get_config()


class ClusterCells:
    def __init__(self, msg):
        self.task_def = msg["body"]
        self.colors = COLOR_POOL.copy()

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
        request = {
            "type": self.task_def["type"],
            "config": self.task_def.get("config", {"resolution": 0.5}),
        }
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getClusters",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        resR = r.json()
        #
        # This is a questionable bit of code, but basically it was a simple way of adjusting the results to the shape expected by the UI
        # Doing this allowed me to use the format function as is.
        # It shouldn't be too taxing, at most O(n of cells), which is well within our time complexity because the taxing part will be clustering.
        #
        resR = pd.DataFrame(resR)
        resR.set_index("_row", inplace=True)
        resR["cluster"] = pd.Categorical(resR.cluster)
        return self._format_result(
            resR, self.task_def["cellSetKey"], self.task_def["cellSetName"]
        )
