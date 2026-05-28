import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder

from ..helpers.get_bucketed_heatmap_cells import get_bucketed_heatmap_cells

from ..config import config
from ..result import Result
from ..tasks import Task
from ..helpers.s3 import get_cell_sets


class GeneExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = self.task_def
        cell_ids = request.get("downsampleSettings", {}).get("cellIds")
        downsample_settings = request.get("downsampleSettings")

        if cell_ids is None and downsample_settings:
            cell_sets = get_cell_sets(self.experiment_id)
            cell_ids = get_bucketed_heatmap_cells(
                downsample_settings["selectedCellSet"],
                downsample_settings["groupedTracks"],
                cell_sets
            )

        request["cellIds"] = cell_ids
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

        response.raise_for_status()
        
        # Response body is already the data as bytes (JSON encoded)
        result = response.content
        return self._format_result(result)
