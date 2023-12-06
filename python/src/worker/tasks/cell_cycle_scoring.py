import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from exceptions import raise_if_error

from ..config import config
from ..result import Result
from ..tasks import Task


class CellCycleScoring(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID
        self.request = msg

    def _format_result(self, result):
        return Result(result, cacheable=False)

    def _format_request(self):
        return { 
            "apiUrl" : config.API_URL, 
            "authJwt" : self.request["Authorization"],
            "experimentId": self.experiment_id
        }


    @xray_recorder.capture("CellCycleScoring.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/CellCycleScoring",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)

