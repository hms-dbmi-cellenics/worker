import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.r_worker_exception import RWorkerException
from ..result import Result
from ..tasks import Task


class GeneExpression(Task):
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
        # raise an exception if an HTTPError if one occurred because otherwise response.json() will fail
        response.raise_for_status()
        result = response.json()

        error = result.get("error", False)
        if error:
            err_message = error.get("message", "")
            err_code = error.get("code", "")
            raise RWorkerException(message=err_message, code=err_code)

        data = result.get("data")
        truncatedExpression = data["truncatedExpression"]
        rawExpression = data["rawExpression"]
        result = {}

        for gene in rawExpression.keys():

            view = rawExpression[gene]
            # can't do summary stats on list with None's
            # casting to np array replaces None with np.nan
            viewnp = np.array(view, dtype=np.float)
            # This is not necessary and is also costly, but I leave it commented as a reminder
            # that this object has integer zeros and floating point for n!=0.
            # expression = [float(item) for item in view]
            mean = float(np.nanmean(viewnp))
            stdev = float(np.nanstd(viewnp))
            result[gene] = {"truncatedExpression": {}, "rawExpression": {}}
            result[gene]["rawExpression"] = {
                "mean": mean,
                "stdev": stdev,
                "expression": view,
            }

            viewTr = truncatedExpression[gene]
            viewnpTr = np.array(viewTr, dtype=np.float)
            minimum = float(np.nanmin(viewnpTr))
            maximum = float(np.nanmax(viewnpTr))
            result[gene]["truncatedExpression"] = {
                "min": minimum,
                "max": maximum,
                "expression": viewTr,
            }

        return self._format_result(result)
