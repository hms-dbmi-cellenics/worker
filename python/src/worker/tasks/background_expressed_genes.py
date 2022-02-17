import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.get_diff_expr_cellsets import get_diff_expr_cellsets
from ..helpers.remove_regex import remove_regex
from ..helpers.s3 import get_cell_sets
from ..result import Result
from . import Task


class GetBackgroundExpressedGenes(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID
        self.pagination = {}

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        # get cell sets from database
        all_cell_sets = get_cell_sets(self.experiment_id)

        first_cell_set_name = self.task_def["cellSet"]
        second_cell_set_name = self.task_def["compareWith"]
        basis_name = self.task_def["basis"]

        baseCells, backgroundCells = get_diff_expr_cellsets(
            basis_name, first_cell_set_name, second_cell_set_name, all_cell_sets
        )

        request = {
            "baseCells": [int(x) for x in baseCells],
            "backgroundCells": [int(x) for x in backgroundCells],
        }

        return request

    @xray_recorder.capture("getBackgroundExpressedGenes.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):

        request = self._format_request()

        # send request to r worker
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getBackgroundExpressedGenes",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise response.json()
        #  will fail
        response.raise_for_status()
        result = response.json()

        error = result.get("error", False)
        if error:
            raise Exception(error)

        return self._format_result({"genes": result["genes"]})
