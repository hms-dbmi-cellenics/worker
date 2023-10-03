import json
import base64

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..helpers.cell_sets_dict import get_cell_sets_dict, subset_cell_sets_dict
from ..result import Result
from ..tasks import Task

class GetNormalizedExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        return Result(result)

    def _format_request(self):
        subset_by = self.task_def["subsetBy"]

        cell_sets = get_cell_sets(self.experiment_id)

        # categories should be ["sample", "louvain", "metadata", "scratchpad"] (in no particular order)
        categories = list(subset_by.keys())

        cell_ids_to_intersect = []

        # There's no subsetting to be done just send a None
        if (all(len(subset_by[category]) == 0 for category in categories)):
            return { "subsetBy": None, "applySubset": False }

        cell_sets_dict = get_cell_sets_dict(cell_sets)

        # Get the sets of cell ids to subset by for each category
        for category in categories:
            if (len(subset_by[category]) == 0):
                continue

            cell_ids = subset_cell_sets_dict(subset_by[category], cell_sets_dict)
            cell_ids_to_intersect.append(cell_ids)

        # Intersect all sets of cell ids from different categories
        cell_ids = cell_ids_to_intersect[0].intersection(*cell_ids_to_intersect[1:])

        return { "subsetBy": list(cell_ids), "applySubset": True }

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

        print("responseDebug")
        print(response)

        print("responseREasonDebug")
        print(response.reason)


        response.raise_for_status()
        result = response.text
        # raise_if_error(result)

        print("resultDebug")
        print(result)

        encoded_chunks = json.loads(result)
        # encoded_chunks = result.get("data")
        
        

        # Decode the base64 chunks and concatenate them
        decoded_data = ''.join([base64.b64decode(chunk).decode('utf-8') for chunk in encoded_chunks])

        return self._format_result(decoded_data)

