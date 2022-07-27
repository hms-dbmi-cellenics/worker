import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..result import Result
from ..tasks import Task


class GetImgPlot(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.task_etag = msg["ETag"]

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result, upload=False)

    def _format_request(self):
        request = {
            "plotType": self.task_def["type"],
            "genes": self.task_def["features"],
            "etag": self.task_etag,
            "plotSubType": self.task_def["plotSubType"],
        }
        return request

    @xray_recorder.capture("GetImgPlot.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getImgPlot",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)
        data = result.get("data")
        # obj = json.loads(data)
        return self._format_result(data)
