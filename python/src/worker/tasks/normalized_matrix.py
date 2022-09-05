import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
import pandas as pd
from exceptions import raise_if_error

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class GetNormalizedExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # convert result to dataframe
        result = pd.DataFrame(result, index=result["_row"], columns = result.keys() - {'_row'})
        
        return Result(result)

    def _format_request(self):
        filterBy = self.task_def["filterBy"]
        applyFilter = all(group.lower() != "all" for group in filterBy["group"])

        if not applyFilter:
            return {
            "filterBy": [],
            "applyFilter": applyFilter,
        }

        # Getting cell ids to subset the seurat object with a group of cells
        cellSets = get_cell_sets(self.experiment_id) 
        filterByCellSets = []
        for cellSet in cellSets:
                if cellSet["key"] in filterBy["group"]:
                    filterByCellSets.append(cellSet)

        filterByCellIds = []
        for cellSet in filterByCellSets:
            for child in cellSet["children"]:
                if child["key"] in filterBy["key"]:
                    filterByCellIds.extend(child["cellIds"])

        return {
            "filterBy": list(set(filterByCellIds)),
            "applyFilter": applyFilter,
        }


    @xray_recorder.capture("GetNormalizedExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/GetNormalizedExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)

