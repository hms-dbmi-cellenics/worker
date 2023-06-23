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

class MarkerHeatmap(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = {"nGenes": self.task_def["nGenes"]}

        downsample_settings = self.task_def["downsampleSettings"]

        selected_cell_set = downsample_settings["selectedCellSet"]
        grouped_tracks = downsample_settings["groupedTracks"]
        selected_points = downsample_settings["selectedPoints"]
        hidden_cell_set_keys = downsample_settings["hiddenCellSets"]

        cell_sets = get_cell_sets(self.experiment_id)

        # There is no max_cells being sent right now, but this allows
        #  us to control this number in the future if we want to
        max_cells = self.task_def.get("maxCells", 1000)

        cell_order = get_heatmap_cell_order(
            selected_cell_set,
            grouped_tracks,
            selected_points,
            hidden_cell_set_keys,
            max_cells,
            cell_sets
        )

        selected_cell_sets = next(set for set in cell_sets if set["key"] == selected_cell_set)
        
        request["cellSets"] = selected_cell_sets
        request["cellIds"] = cell_order
        return request, cell_order

    @xray_recorder.capture("MarkerHeatmap.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request, cell_order = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/runMarkerHeatmap",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        json_response = response.json()
        raise_if_error(json_response)
        result = json_response.get("data")

        result["cellOrder"] = cell_order

        return self._format_result(result)
