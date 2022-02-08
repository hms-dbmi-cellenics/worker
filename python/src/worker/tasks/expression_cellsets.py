import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..result import Result
from ..tasks import Task


class GetExpressionCellSets(Task):
    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _construct_request(self):
        request = self.task_def
        request["config"]["experimentId"] = config.EXPERIMENT_ID
        request["config"]["apiUrl"] = config.API_URL
        return request


    @xray_recorder.capture("getExpressionCellSets.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._construct_request()

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getExpressionCellSets",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise
        #  r.json() will fail
        r.raise_for_status()
        result = r.json()
        return self._format_result(result)
