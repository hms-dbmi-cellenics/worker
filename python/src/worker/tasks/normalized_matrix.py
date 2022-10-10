import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
import pandas as pd
from exceptions import raise_if_error

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task

# If keys_to_subset
def get_cell_ids(category, cell_sets, keys_to_subset = []):
    # Get all cell sets in the category
    cell_sets_in_category = next(cell_class for cell_class in cell_sets if cell_class["key"] == category)["children"]

    cell_ids = set()

    for cell_set in cell_sets_in_category:
        if (keys_to_subset == True or cell_set["key"] in keys_to_subset):
            cell_ids = cell_ids.union(set(cell_set["cellIds"]))


class GetNormalizedExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # convert result to dataframe
        result = pd.DataFrame(result, index=result["_row"], columns = result.keys() - {'_row'})
        
        return Result(result)

    def _format_request(self):
        subset_by = self.task_def["subsetBy"]
        # apply_filter = all(group.lower() != "all" for group in subset_by["group"])

        # sample: [],
        # louvain: [],
        # metadata: [],
        # scratchpad: [],

        cell_sets = get_cell_sets(self.experiment_id)

        categories = ["sample", "louvain", "metadata", "scratchpad"]

        cell_ids_to_intersect = []

        if (all(len(subset_by[category]) == 0 for category in categories)):
            pass
            # TODO: return all_cells

        for category in categories:
            if (len(subset_by[category])):
                continue

            cell_ids = get_cell_ids(subset_by[category], category, cell_sets)

            cell_ids_to_intersect.append(cell_ids)

        cell_ids = cell_ids_to_intersect[0].intersection(*cell_ids_to_intersect[1:])

        return {
            "subsetBy": list(cell_ids)
        }
        # if not apply_filter:
        #     return {
        #     "filterBy": [],
        #     "applyFilter": apply_filter,
        # }

        # # Getting cell ids to subset the seurat object with a group of cells
        # cellSets = get_cell_sets(self.experiment_id) 
        # filterByCellSets = []
        # for cellSet in cellSets:
        #     if cellSet["key"] in subset_by["group"]:
        #         filterByCellSets.append(cellSet)

        # filterByCellIds = []
        # for cellSet in filterByCellSets:
        #     for child in cellSet["children"]:
        #         if child["key"] in subset_by["key"]:
        #             filterByCellIds.extend(child["cellIds"])

        # return {
        #     "filterBy": list(set(filterByCellIds)),
        #     "applyFilter": apply_filter,
        # }


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

