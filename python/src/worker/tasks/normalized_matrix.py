import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
import pandas as pd
from exceptions import raise_if_error
from io import StringIO
from collections import defaultdict

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task

# Move all cell_sets data into a dict with { cellSetKey: cellSetIds } for faster access
def get_cell_sets_dict(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        for cell_set in cell_class["children"]:
            cell_sets_dict[cell_set["key"]] = cell_set["cellIds"]

    return cell_sets_dict

# Get all cell sets that match the subset_keys
def get_cell_ids(subset_keys, cell_sets_dict):
    cell_ids = set()

    for subset_key in subset_keys:
        cell_ids = cell_ids.union(set(cell_sets_dict[subset_key]))
    
    return cell_ids


class GetNormalizedExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # convert result to dataframe
        print("Reading csv into pandas dataframe")
        csvStringIO = StringIO(result)
        result = pd.read_csv(csvStringIO, sep=",", header=None, dtype=defaultdict(lambda: int))

        print("Finished creating pandas dataframe")
        # result = pd.DataFrame(result, index=result["_row"], columns = result.keys() - {'_row'})
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
            
            cell_ids = get_cell_ids(subset_by[category], cell_sets_dict)
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

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        return self._format_result(data)

