import json
from config import get_config
from result import Result
import pandas
import requests

from helpers.dynamo import get_item_from_dynamo
from helpers.find_cells_by_set_id import find_cells_by_set_id
from helpers.find_cell_ids_in_same_hierarchy import find_cell_ids_in_same_hierarchy

config = get_config()


class DifferentialExpression:
    def __init__(self, msg, adata):
        self.adata = adata
        self.task_def = msg["body"]
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        result = result.to_dict(orient="records")

        # JSONify result.
        result = json.dumps({"rows": result})

        # Return a list of formatted results.
        return [Result(result)]

    # Marks with a special string the cells that are not included in the experiment
    def mark_cells_not_in_basis_set(self, raw_adata, resp):
        basis_cell_set_name = self.task_def.get("basis")

        # Ignore basis if it is all cells in the experiment or not specified (default behaviour for previous versions)
        if not basis_cell_set_name or ("all" in basis_cell_set_name.lower()):
            return set()

        basis_cell_set = set(find_cells_by_set_id(basis_cell_set_name, resp))

        # Mark all cells that are not included in the sample as FilteredOut, these will be ignored in the rest of the experiment
        raw_adata.obs["condition"].loc[
            ~raw_adata.obs["cell_ids"].isin(basis_cell_set)
        ] = "FilteredOut"

    # Get cells values for the cell set.
    def get_cells_in_set(self, name, resp, first_cell_set_name):
        cells = []

        # If "rest", then get all cells in the same hierarchy as the first cell set that arent part of "first"
        if "rest" in name.lower():
            print("found rest")
            cells = find_cell_ids_in_same_hierarchy(first_cell_set_name, resp)
        else:
            cells = find_cells_by_set_id(name, resp)

        return cells

    def mark_cells_of_set_in_data(self, label, cells, raw_adata):
        # Append the current label to the already existing one in "condition",
        # this way when getting the cells for one cell set
        # we can ignore those that are intersected with the other one
        # since these would be problematic for the factor() function in the r worker.
        raw_adata.obs["condition"].loc[
            (raw_adata.obs["cell_ids"].isin(cells))
            & (~raw_adata.obs["condition"].eq("FilteredOut"))
        ] += label

    def compute(self):
        # get the top x number of genes to load:
        n_genes = self.task_def.get("maxNum", None)

        # get cell sets from database
        resp = get_item_from_dynamo(self.experiment_id, "cellSets")

        # use raw values for this task
        raw_adata = self.adata.raw.to_adata()

        # create a series to hold the conditions
        raw_adata.obs["condition"] = ""

        # mark all cells that aren't in the basis sample to be filtered out
        self.mark_cells_not_in_basis_set(raw_adata, resp)

        first_cell_set_name = self.task_def["cellSet"]
        second_cell_set_name = self.task_def["compareWith"]

        # mark cells of first set
        first_cell_set = find_cells_by_set_id(first_cell_set_name, resp)
        self.mark_cells_of_set_in_data("first", first_cell_set, raw_adata)

        # mark cells of second set
        if (
            second_cell_set_name == "background"
            or "all" in second_cell_set_name.lower()
        ):
            raw_adata.obs["condition"].loc[raw_adata.obs["condition"] == ""] = "second"

        else:
            second_cell_set = self.get_cells_in_set(
                second_cell_set_name, resp, first_cell_set_name
            )

            self.mark_cells_of_set_in_data("second", second_cell_set, raw_adata)

        # create request from the corresponding marked cells
        request = {
            "baseCells": raw_adata.obs.index[
                raw_adata.obs["condition"] == "first"
            ].tolist(),
            "backgroundCells": raw_adata.obs.index[
                raw_adata.obs["condition"] == "second"
            ].tolist(),
        }

        # send request to r worker
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        result = pandas.DataFrame.from_dict(r.json())

        result.dropna(inplace=True)

        # get top x most significant results, if parameter was supplied
        if n_genes:
            result = result.nsmallest(n_genes, ["abszscore"])

        return self._format_result(result)
