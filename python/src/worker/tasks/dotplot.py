import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.r_worker_exception import RWorkerException
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class DotPlot(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):

        # getting cell ids for the groups we want to display.
        cellSets = get_cell_sets(self.experiment_id)

        # Getting the cell ids for subsetting the seurat object with a group of cells.
        groupByCellSet = [
            cellSet
            for cellSet in cellSets
            if cellSet["key"] == self.task_def["groupBy"]
        ][0]

        filterBy = self.task_def["filterBy"]
        applyFilter = filterBy["group"].lower() != "all"
        filterByCellSet = groupByCellSet

        if applyFilter:
            children = [
                cellSet for cellSet in cellSets if cellSet["key"] == filterBy["group"]
            ][0]["children"]
            filterByCellSet = [
                child for child in children if child["key"] == filterBy["key"]
            ][0]

        request = {
            "useMarkerGenes": self.task_def["useMarkerGenes"],
            "numberOfMarkers": self.task_def["numberOfMarkers"],
            "customGenesList": self.task_def["customGenesList"],
            "groupBy": groupByCellSet,
            "filterBy": filterByCellSet,
            "applyFilter": applyFilter,
        }

        return request

    @xray_recorder.capture("DotPlot.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):

        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/runDotPlot",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        # raise an exception if an HTTPError if one occurred because otherwise response.json() will fail
        response.raise_for_status()
        result = response.json()

        error = result.get("error", False)
        if error:
            user_msg = result.get("user_msg", "")
            raise RWorkerException(user_msg=user_msg, error=error)

        return self._format_result(result)
