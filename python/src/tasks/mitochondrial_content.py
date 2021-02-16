import json
from config import get_config
from result import Result
import requests

config = get_config()


class GetMitochondrialContent:
    def __init__(self, msg):
        self.task_def = msg["body"]

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # the cells to get MT-content data for
        cells = self.task_def["cells"]

        request = {"cells": cells}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        
        resultR = r.json()
        result = {}
        for i in range(len(resultR['cells_id'])):
            cell = resultR['cells_id'][i]
            doublet_scores = resultR['percent.mt'][i]

            result[cell] = {
                "percent-MT": doublet_scores,
            }

        return self._format_result(result)
