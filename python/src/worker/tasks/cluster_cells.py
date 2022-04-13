import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.color_pool import COLOR_POOL
from ..result import Result
from ..tasks import Task


class ClusterCells(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.colors = COLOR_POOL.copy()
        self.request = msg

    def _format_result(self, result):
        return Result(result, cacheable=False)

    def _format_request(self):
        resolution = self.task_def["config"].get("resolution", 0.5)

    # add apiUrl and authJwt to req to be able to patch in R worker
        request = {
            "type": self.task_def["type"],
            "config": {"resolution" : resolution},
            "apiUrl" : config.API_URL,
            "authJwt" : self.request["Authorization"]
        }
        
        print(request)
        
        return request

    @xray_recorder.capture("ClusterCells.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):

        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getClusters",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)
