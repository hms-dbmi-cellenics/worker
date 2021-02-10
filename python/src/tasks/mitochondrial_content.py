import json
from config import get_config
from result import Result
import requests

config = get_config()


class MitochondrialContent:
    def __init__(self, msg, adata):
        self.adata = adata
        self.task_def = msg["body"]

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # the genes to get expression data for
        cells = self.task_def["cells"]

        request = {"cells": cells}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        
        resultR = r.json()
        result = {}
        for i in range(len(resultR['_row'])):
            cell = resultR['_row'][i]
            doublet_scores = resultR['percent.mt'][i]

            result[cell] = {
                "percent-MT": doublet_scores,
            }

        return self._format_result(result)
