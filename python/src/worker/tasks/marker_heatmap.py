import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.process_gene_expression import process_gene_expression
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
        raise_if_error(json_response)

        data = json_response.get("data")
        result = {
            "data": process_gene_expression(data),
            "order": list(data["rawExpression"]),
        }

        return self._format_result(result)
