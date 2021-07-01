import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..result import Result
from ..tasks import Task


class GetDoubletScore(Task):
    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    @xray_recorder.capture("DoubletScore.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):

        # Retrieve the MitochondrialContent of all the cells
        request = {}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getDoubletScore",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        # The values are ordered by cells id
        # The result contains a list with the doublet scores values
        result = r.json()
        return self._format_result(result)
