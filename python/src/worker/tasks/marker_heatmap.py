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

def get_heatmap_cell_order(n_genes, cell_set_key, grouped_tracks, selected_points, hidden_cell_set_keys, cell_sets):


class MarkerHeatmap(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result)

    def _format_request(self):
        request = {"nGenes": self.task_def["nGenes"], "cellIds": self.task_def["cellIds"]}

        cellSetKey = self.task_def["cellSetKey"]

        cell_sets = get_cell_sets(self.experiment_id)

        n_genes = self.task_def["nGenes"]
        cell_set_key = self.task_def["cellSetKey"]
        grouped_tracks = self.task_def["groupByClasses"]
        selected_points = self.task_def["selectedPoints"]
        hidden_cell_set_keys = self.task_def["hiddenCellSetKeys"]

        cell_order = get_heatmap_cell_order(
            n_genes,
            cell_set_key,
            grouped_tracks,
            selected_points,
            hidden_cell_set_keys,
            cell_sets
        )

        for set in cell_sets:
            if set["key"] == cellSetKey:
                cell_sets = set
                break

        request["cellSets"] = cell_sets
        request["cellIds"] = cell_order
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

        response.raise_for_status()
        json_response = response.json()
        raise_if_error(json_response)
        result = json_response.get("data")
        return self._format_result(result)
