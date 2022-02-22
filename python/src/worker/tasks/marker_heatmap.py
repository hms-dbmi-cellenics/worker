import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.worker_exception import WorkerException
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class MarkerHeatmap(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = {"nGenes": self.task_def["nGenes"]}

        cellSetKey = self.task_def["cellSetKey"]

        cellSets = get_cell_sets(self.experiment_id)

        for set in cellSets:
            if set["key"] == cellSetKey:
                cellSets = set
                break

        request["cellSets"] = cellSets
        return request

    @xray_recorder.capture("MarkerHeatmap.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/runMarkerHeatmap",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError occurred
        # as otherwise response.json() will fail
        response.raise_for_status()
        json_response = response.json()

        error = json_response.get("error", False)
        if error:
            user_message = error.get("user_message", "")
            err_code = error.get("error_code", "")
            raise WorkerException(err_code, user_message)

        response_data = json_response.get("data")

        truncatedExpression = response_data["truncatedExpression"]
        rawExpression = response_data["rawExpression"]
        result = {}
        data = {}
        order = []

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
            data[gene] = {"truncatedExpression": {}, "rawExpression": {}}
            data[gene]["rawExpression"] = {
                "mean": mean,
                "stdev": stdev,
                "expression": view,
            }

            viewTr = truncatedExpression[gene]
            viewnpTr = np.array(viewTr, dtype=np.float)
            minimum = float(np.nanmin(viewnpTr))
            maximum = float(np.nanmax(viewnpTr))
            data[gene]["truncatedExpression"] = {
                "min": minimum,
                "max": maximum,
                "expression": viewTr,
            }
            order.append(gene)

        result["data"] = data
        result["order"] = order
        return self._format_result(result)
