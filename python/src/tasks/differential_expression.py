import json
from config import get_config
from result import Result
import pandas
import requests

from helpers.dynamo import get_item_from_dynamo
from helpers.find_cells_by_set_id import find_cells_by_set_id

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

    def compute(self):
        # the cell set to compute differential expression on
        cell_sets = {
            'first': self.task_def["cellSet"],
            'second': self.task_def["compareWith"]
        }
        
        # get the top x number of genes to load:
        n_genes = self.task_def.get("maxNum", None)

        # get cell sets from database
        resp = get_item_from_dynamo(self.experiment_id, "cellSets")

        # use raw values for this task
        raw_adata = self.adata.raw.to_adata()

        # create a series to hold the conditions
        raw_adata.obs["condition"] = None

        # fill in values appropriately. if `rest`, fill in all NaN
        # values with the label, as the other one will be a cell set
        # and override the appropriate values. if between two cell sets,
        # make sure the cells not in either will be marked with `other`.
        for label, name in cell_sets.items():
            if name == "rest" or "all" in name:
                raw_adata.obs["condition"].fillna(value=label, inplace=True)
            else:
                cells = find_cells_by_set_id(name, resp)

                raw_adata.obs["condition"].loc[
                    raw_adata.obs["cell_ids"].isin(cells)
                ] = label
        
        raw_adata.obs["condition"].fillna(value="other", inplace=True)

        raw_adata.obs = raw_adata.obs.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)

        request = {
            "baseCells": raw_adata.obs.index[
                raw_adata.obs["condition"] == "first"
            ].tolist(),
            "backgroundCells": raw_adata.obs.index[
                raw_adata.obs["condition"] == "second"
            ].tolist(),
        }

        print(request)

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
