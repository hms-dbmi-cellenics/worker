import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..result import Result
from ..tasks import Task


class GetMitochondrialContent(Task):
    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    @xray_recorder.capture("GetMitochondrialContent.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def _format_request(self):
        return {}

    def compute(self):
        # Retrieve the MitochondrialContent of all the cells
        request = self._format_request()
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError occurred
        # as otherwise response.json() will fail
        response.raise_for_status()
        # The values are ordered by cells id
        # The result contains a list with the MT-content values
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)
