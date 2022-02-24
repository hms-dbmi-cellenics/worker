import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.get_diff_expr_cellsets import get_diff_expr_cellsets
from ..helpers.remove_regex import remove_regex
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class DifferentialExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID
        self.pagination = {}

        if "pagination" in msg:
            self.pagination = msg["pagination"]

    def _format_result(self, result):
        # Return a list of formatted results.

        return Result({"total": result["full_count"], "rows": result["gene_results"]})

    def _format_request(self):
        # get cell sets from database
        cell_sets = get_cell_sets(self.experiment_id)
        first_cell_set_name = self.task_def["cellSet"]
        second_cell_set_name = self.task_def["compareWith"]
        basis_name = self.task_def["basis"]

        baseCells, backgroundCells = get_diff_expr_cellsets(
            basis_name, first_cell_set_name, second_cell_set_name, cell_sets
        )

        request = {
            "baseCells": [int(x) for x in baseCells],
            "backgroundCells": [int(x) for x in backgroundCells],
            "genesOnly": self.task_def.get("genesOnly", False),
            "comparisonType": self.task_def.get("comparisonType", "within"),
        }

        if self.pagination:
            request["pagination"] = self.pagination

        if (
            "filters" in self.pagination
            and self.pagination["filters"][0]["type"] == "text"
        ):
            gene_filter = self.pagination["filters"][0]["expression"]
            request["geneNamesFilter"] = remove_regex(gene_filter)

        return request

    @xray_recorder.capture("DifferentialExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):

        request = self._format_request()

        # send request to r worker
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)
