import json

import backoff
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

        # getting cell ids for the groups we want to display.
        cellSets = get_cell_sets(self.experiment_id)
        setNames = [set["key"] for set in cellSets]
    
        # Getting the cell ids for subsetting the seurat object with a group of cells.
        cellSets = cellSets[setNames.index(self.task_def["groupBy"])]
        filterBy = self.task_def["filterBy"]
        filterByAll = filterBy['group'].lower() == "all"

        if filterByAll:
            filterBy = cellSets
        else:
            children = cellSets[setNames.index(filterBy['group'])]["children"]
            childrenNames = [set["key"] for set in children]
            filterBy = childrenNames[childrenNames.index(filterBy['key'])]

        request = {
            "useMarkerGenes": self.task_def["useMarkerGenes"],
            "numberOfMarkers": self.task_def["numberOfMarkers"],
            "customGenesList": self.task_def["customGenesList"],
            "cellSets": cellSets,
            "filterBy": filterBy,
            "filterByAll": filterByAll,
        }

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/runDotPlot",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        result = r.json()

        return self._format_result(result)
