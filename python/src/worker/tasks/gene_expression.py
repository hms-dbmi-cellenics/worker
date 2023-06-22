import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.process_gene_expression import process_gene_expression
from ..helpers.get_heatmap_cell_order import get_heatmap_cell_order
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class GeneExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = self.task_def

        # If downsampled then calculate cell order from settings
        # And add it to the request
        if (self.task_def["downsampled"] == True):
            cell_sets = get_cell_sets(self.experiment_id)

            downsample_settings = self.task_def["downsampleSettings"]

            cell_set_key = downsample_settings["selectedCellSet"]
            grouped_tracks = downsample_settings["groupedTracks"]
            selected_points = downsample_settings["selectedPoints"]
            hidden_cell_set_keys = downsample_settings["hiddenCellSets"]
            # cell_set_key = downsample_settings["cellSetKey"]
            # grouped_tracks = downsample_settings["groupByClasses"]
            # selected_points = downsample_settings["selectedPoints"]
            # hidden_cell_set_keys = downsample_settings["hiddenCellSetKeys"]

            # There is no max_cells being sent right now, but this allows
            #  us to control this number in the future if we want to
            max_cells = downsample_settings.get("maxCells", 1000)

            cell_order = get_heatmap_cell_order(
                cell_set_key,
                grouped_tracks,
                selected_points,
                hidden_cell_set_keys,
                max_cells,
                cell_sets
            )
            
            request["cellIds"] = cell_order
        return request

    @xray_recorder.capture("GeneExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()
        
        cell_order = request.get("cellIds")

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/runExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)
        result = result.get("data")

        if cell_order:
            result["cellOrder"] = cell_order

        return self._format_result(result)
