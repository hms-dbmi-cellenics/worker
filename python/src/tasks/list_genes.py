import pandas as pd
from result import Result
import requests
from config import get_config
import json

config = get_config()


class ListGenes:
    def __init__(self, msg, adata):
        self.task_def = msg["body"]
        self.adata = adata

    def _format_result(self, result, total):
        # convert result to list of row dicts
        result = result.to_dict(orient="records")

        # JSONify result.
        result = json.dumps({"total": total, "rows": result})

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        request = self.task_def
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/listGenes",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        resR = r.json()
        # Convert to dataframe to prepare the data for the UI
        resR = pd.DataFrame(resR)
        total = 0
        if len(resR) > 0:
            total = resR["full_count"][0]
        resR = resR.drop("full_count", axis=1)
        # total returns numpy int64, convert to integer
        # for serialization to JSON
        total = int(total)
        return self._format_result(resR, total=total)