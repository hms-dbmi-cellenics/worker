import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.r_worker_exception import RWorkerException
from ..result import Result
from ..tasks import Task


class GetExpressionCellSets(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.request = msg

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _construct_request(self):
        request = self.task_def

        request["config"] = {
            "experimentId": config.EXPERIMENT_ID,
            "apiUrl": config.API_URL,
            "authJwt": self.request["Authorization"],
        }
        return request

    @xray_recorder.capture("getExpressionCellSet.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._construct_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getExpressionCellSet",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise
        #  response.json() will fail
        response.raise_for_status()
        result = response.json()

        error = result.get("error", False)
        if error:
            err_message = error.get("message", "")
            err_code = error.get("code", "")
            raise RWorkerException(message=err_message, code=err_code)

        data = result.get("data")

        return self._format_result(data)
