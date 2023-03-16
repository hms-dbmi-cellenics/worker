import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
# from ..helpers.cell_sets_dict import get_cell_sets_dict
from ..result import Result
from ..tasks import Task
from ..helpers.s3 import get_embedding, get_cell_sets

import os

# Move all cell_sets data into a dict
# TODO: MOVE THIS INTO HELPERS
def get_cell_sets_dict_sctype(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        if cell_class["name"].startswith("ScType-"):
            cell_sets_dict[cell_class["name"]] = cell_class["children"]
        else:
            cell_sets_dict[cell_class["key"]] = cell_class["children"]

    return cell_sets_dict

class DownloadAnnotSeuratObject(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        return Result(result)

    def _format_request(self):

        cell_sets = get_cell_sets(self.experiment_id)
        cell_sets_dict = get_cell_sets_dict_sctype(cell_sets)

        embedding_etag = self.task_def["embedding"]["ETag"]
        embedding = get_embedding(embedding_etag, format_for_r=True)

        request = {
            "expId": config.EXPERIMENT_ID,
            "bucket": config.RESULTS_BUCKET,
            "embedding": embedding,
            "embedding_settings": {
                "method": self.task_def["embedding"]["method"],
                "methodSettings": self.task_def["embedding"]["methodSettings"]
            },
            "cellSets": cell_sets_dict,
        }

        return request

    @xray_recorder.capture("DownloadAnnotSeuratObject.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/DownloadAnnotSeuratObject",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        # path to the r.rds file
        rds_file_path = "/R/r.rds"

        # check if the file exists
        if os.path.exists(rds_file_path):
            print("File r.rds exists")
        else:
            print("File not found.")

        return self._format_result(rds_file_path)
    