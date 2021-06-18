import json
import backoff
from config import config
from result import Result
from tasks import Task
import requests
from aws_xray_sdk.core import xray_recorder


class GetMitochondrialContent(Task):
    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    @xray_recorder.capture('GetMitochondrialContent.compute')
    @backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=30)
    def compute(self):
        # Retrieve the MitochondrialContent of all the cells
        request = {}
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/getMitochondrialContent",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        # The values are ordered by cells id
        # The result contains a list with the MT-content values
        result = r.json()
        return self._format_result(result)
