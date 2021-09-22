import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class DotPlot(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)
        # Return a list of formatted results.
        return [Result(result)]

    @xray_recorder.capture("DotPlot.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = {"nGenes": self.task_def["nGenes"], "type":self.task_def["type"], "config":self.task_def["config"]}
        
        if self.task_def["type"] is "custom":
            request["genes"] = self.task_def["genes"]
        else:
            request["nGenes"] = self.task_def["nGenes"]

        typeOfSets = self.task_def["typeOfSets"]

        cellSets = get_cell_sets(self.experiment_id)

        setNames = [set["key"] for set in cellSets]
        request["cellSets"] = cellSets[setNames.index(typeOfSets)]

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/runDotPlot",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        result = r.json()

        return self._format_result(result)
