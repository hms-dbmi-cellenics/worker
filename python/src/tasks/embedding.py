import json
import requests
import backoff
from config import get_config
from result import Result
from aws_xray_sdk.core import xray_recorder

config = get_config()


class ComputeEmbedding:
    def __init__(self, msg):
        self.task_def = msg["body"]

    def _format_result(self, raw):
        # JSONify result.
        raw_result = json.dumps(raw)

        # Return a list of formatted results.
        return [Result(raw_result), Result(raw_result)]

    @xray_recorder.capture('ComputeEmbedding.compute')
    @backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=30)
    def compute(self):
        request = {"type": self.task_def["type"], "config": self.task_def["config"]}

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getEmbedding",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # The index order relies on cells_id in an ascending form. The order is made in the R part. 
        result = r.json()
        return self._format_result(result)
