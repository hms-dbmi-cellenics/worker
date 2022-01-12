import json

import backoff
import pandas
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.find_cell_ids_in_same_hierarchy import (
    find_all_cell_ids_in_cell_sets, find_cell_ids_in_same_hierarchy)
from ..helpers.find_cells_by_set_id import find_cells_by_set_id
from ..helpers.remove_regex import remove_regex
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class DifferentialExpression(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID
        self.pagination = {}
        if "pagination" in msg:
            self.pagination = msg["pagination"]

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result({"total": result["full_count"], "rows": result["gene_results"]})

    def _format_request(self):
        # get cell sets from database
        resp = get_cell_sets(self.experiment_id)

        first_cell_set_name = self.task_def["cellSet"]
        second_cell_set_name = self.task_def["compareWith"]
        basis = self.task_def["basis"]

        # Check if the comparsion is between all the cells or within a cluster
        if not basis or ("all" in basis.lower()):
            # In case that we are not going to filter by a cluster we remain the set
            #  empty
            filtered_set = set()
        else:
            # In the case that we filter by a cluster we keep the cells id in a set
            filtered_set = set(find_cells_by_set_id(basis, resp))

        # mark cells of first set
        first_cell_set = find_cells_by_set_id(first_cell_set_name, resp)

        # mark cells of second set
        # check if the second set is composed by the "All other cells"
        if (
            second_cell_set_name == "background"
            or "all" in second_cell_set_name.lower()
        ):
            # Retrieve all the cells (not necessary at the same hierachy level)
            complete_cell_set = set(find_all_cell_ids_in_cell_sets(resp))
            # Filter with those that are not in the first cell set
            second_cell_set = [
                item for item in complete_cell_set if item not in first_cell_set
            ]
        else:
            # In the case that we compare with specific cell set, we just look for
            #  the cell directly
            second_cell_set = self.get_cells_in_set(
                second_cell_set_name, resp, first_cell_set_name
            )
            # Check any possible intersect cells
            inter_cell_set = set(first_cell_set).intersection(set(second_cell_set))
            first_cell_set = [
                item for item in first_cell_set if item not in inter_cell_set
            ]
            second_cell_set = [
                item for item in second_cell_set if item not in inter_cell_set
            ]

        # Keep only the cell_set that are on the specify basis (in the case that we
        # are not in the "All" analysis)
        if len(filtered_set) > 0:
            second_cell_set = [it for it in second_cell_set if it in filtered_set]
            first_cell_set = [it for it in first_cell_set if it in filtered_set]

        # Check if the first cell set is empty
        if len(first_cell_set) == 0:
            raise Exception("No cells id fullfills the 1st cell set.")

        # Check if the second cell set is empty
        if len(second_cell_set) == 0:
            raise Exception("No cells id fullfills the 2nd cell set.")

        request = {
            "baseCells": [int(x) for x in first_cell_set],
            "backgroundCells": [int(x) for x in second_cell_set],
        }

        if self.pagination:
            request["pagination"] = self.pagination

        if "filters" in self.pagination and self.pagination["filters"][0]["type"] == "text":
            gene_filter = self.pagination["filters"][0]["expression"]
            request["geneNamesFilter"] = remove_regex(gene_filter)

        return request

    # Get cells values for the cell set.
    def get_cells_in_set(self, name, resp, first_cell_set_name):
        cells = []

        # If "rest", then get all cells in the same hierarchy as the first cell set
        #  that arent part of "first"
        if "rest" in name.lower():
            cells = find_cell_ids_in_same_hierarchy(first_cell_set_name, resp)
        else:
            cells = find_cells_by_set_id(name, resp)

        return cells

    @xray_recorder.capture("DifferentialExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        
        request = self._format_request()
        
        # send request to r worker
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise r.json()
        #  will fail
        r.raise_for_status()
        r = r.json()
        
        return self._format_result(r)
