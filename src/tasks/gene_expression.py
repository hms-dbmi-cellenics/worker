import json
import numpy as np
from config import get_config
from result import Result

config = get_config()


class GeneExpression:
    def __init__(self, msg, adata):
        self.adata = adata
        self.task_def = msg["body"]
        self.experiment_id = msg["experimentId"]

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # the genes to get expression data for
        genes = self.task_def["genes"]

        # whether to perform feature scaling (defaults to False)
        scale = self.task_def.get("scale", False)

        # compute data on raw matrix
        raw_adata = self.adata.raw.to_adata()
        raw_adata = raw_adata.copy()
        raw_adata.X = raw_adata.X.toarray()

        # create a proper ordering of cells by increasing IDs
        obs_copy = raw_adata.obs.copy()
        obs_copy.sort_values(by=["cell_ids"], inplace=True)
        obs_copy = obs_copy.dropna()
        cell_list = obs_copy.index.tolist()

        # try to find all genes in the list
        raw_adata.var["genes_to_compute"] = raw_adata.var.index.isin(genes)

        # this orders and filters the matrix correctly for both cells and genes
        raw_adata = raw_adata[
            cell_list, raw_adata.var["genes_to_compute"],
        ]

        # if feature scaling is desired, perform that now
        # we do not do this now.
        if scale:
            pass

        result = {}

        for gene in raw_adata.var.index:
            view = raw_adata[:, raw_adata.var.index == gene]

            expression = view.X.flatten().tolist()
            minimum = float(np.amin(view.X))
            maximum = float(np.amax(view.X))

            result[gene] = {
                "min": minimum,
                "max": maximum,
                "expression": expression,
            }

        return self._format_result(result)
