import json
import time

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..result import Result
from ..tasks import Task


class GeneExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = self.task_def
        return request

    @xray_recorder.capture("GeneExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()
        
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/runExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        
        # Response body is already the data as bytes (JSON encoded)
        result = response.content
        return self._format_result(result)
