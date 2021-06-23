import pandas as pd
import json
import requests
import backoff
from result import Result
from helpers.color_pool import COLOR_POOL
from config import get_config
from aws_xray_sdk.core import xray_recorder
from natsort import natsorted

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
            "type": "cellSets",
            "children": [],
        }
        for cluster in natsorted(raw["cluster"].cat.categories):
            view = raw[raw.cluster == cluster]["cell_ids"]
            cell_set["children"].append(
                {
                    "key": f"{cell_set_key}-{cluster}",
                    "name": f"Cluster {cluster}",
                    "rootNode": False,
                    "type": "cellSets",
                    "color": self.colors.pop(0),
                    "cellIds": [int(id_) for id_ in view.tolist()],
                }
            )
        return [Result(json.dumps(cell_set), cacheable=False)]

    @xray_recorder.capture('ClusterCells.compute')
    @backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=30)
    def compute(self):
        resolution = self.task_def["config"].get("resolution",0.5)

        request = {
            "type": self.task_def["type"],
            "config": self.task_def.get("config", {"resolution": resolution}),
        }

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getClusters",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        resR = r.json()

        # This is a questionable bit of code, but basically it was a simple way of adjusting the results to the shape
        # expected by the UI Doing this allowed me to use the format function as is. It shouldn't be too taxing,
        # at most O(n of cells), which is well within our time complexity because the taxing part will be clustering.
        resR = pd.DataFrame(resR)
        resR.set_index("_row", inplace=True)
        resR["cluster"] = pd.Categorical(resR.cluster)
        return self._format_result(
            resR, self.task_def["cellSetKey"], self.task_def["cellSetName"]
        )
