import json
from config import get_config
from result import Result
import numpy as np
import requests

config = get_config()


class GetMitochondrialContent:
    def __init__(self):
        self.task_def = "MT-content"

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # Retrieve the MitochondrialContent of all the cells
        request = {}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # Format result with the same structure than Gene Expression. The values are ordered by the default order
        # of the Embedding in the Seurat object
        # {"MT-content": {'min': , X 'max': , X 'mean': , X 'stdev': , X 'values': [...]}

        result = {}
        values = r.json()
        minimum = float(np.amin(values))
        maximum = float(np.amax(values))
        mean = float(np.mean(values))
        stdev = float(np.std(values))

        result[self.task_def] = {
            "min": minimum,
            "max": maximum,
            "mean": mean,
            "stdev": stdev,
            "values": values,
        }

        return self._format_result(result)

