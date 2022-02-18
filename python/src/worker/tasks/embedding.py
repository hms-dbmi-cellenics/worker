import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.r_worker_exception import RWorkerException
from ..result import Result
from ..tasks import Task


class GetEmbedding(Task):
    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = {
            "type": self.task_def["type"],
            "config": self.task_def["config"],
        }
        return request

    @xray_recorder.capture("ComputeEmbedding.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getEmbedding",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise response.json() will fail
        response.raise_for_status()
        # The index order relies on cells_id in an ascending form. The order is made in the R part.
        result = response.json()

        error = result.get("error", False)
        if error:
            user_message = error.get("user_message", "")
            err_code = error.get("code", "")
            raise RWorkerException(user_message, err_code)

        data = result.get("data")

        return self._format_result(data)
