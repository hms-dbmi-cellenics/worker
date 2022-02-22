import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.worker_exception import WorkerException
from ..result import Result
from ..tasks import Task


class GetDoubletScore(Task):
    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    @xray_recorder.capture("DoubletScore.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def _format_request(self):
        return {}

    def compute(self):

        # Retrieve the MitochondrialContent of all the cells
        request = self._format_request()
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getDoubletScore",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError occurred
        # as otherwise response.json() will fail
        response.raise_for_status()

        # The values are ordered by cells id
        # The result contains a list with the doublet scores values
        result = response.json()

        error = result.get("error", False)
        if error:
            user_message = error.get("user_message", "")
            err_code = error.get("error_code", "")
            raise WorkerException(err_code, user_message)

        data = result.get("data")

        return self._format_result(data)
