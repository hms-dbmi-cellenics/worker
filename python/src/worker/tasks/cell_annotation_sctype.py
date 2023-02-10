import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


# Move all cell_sets data into a dict
def get_cell_sets_dict(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        cell_sets_dict[cell_class["key"]] = cell_class["children"]

    return cell_sets_dict

class ScTypeAnnotate(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID
        self.request = msg

    def _format_result(self, result):
        return Result(result, cacheable=False)

    def _format_request(self):
        # get cell sets from database
        cell_sets = get_cell_sets(self.experiment_id)
        cell_sets_dict = get_cell_sets_dict(cell_sets)

        species = self.task_def["species"]
        tissue = self.task_def["tissue"]

        return { "cellSets": cell_sets_dict, "species": species, "tissue": tissue, "apiUrl" : config.API_URL, "authJwt" : self.request["Authorization"]}


    @xray_recorder.capture("ScTypeAnnotate.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/ScTypeAnnotate",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)

