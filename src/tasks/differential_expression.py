import json
from config import get_config
from result import Result
import diffxpy.api as de

from helpers.dynamo import get_item_from_dynamo
from helpers.find_cells_by_set_id import find_cells_by_set_id

config = get_config()


class DifferentialExpression:
    def __init__(self, msg, adata):
        self.adata = adata
        self.task_def = msg["body"]
        self.experiment_id = msg["experimentId"]

    def _format_result(self, result):
        result = result.to_dict(orient="records")

        # JSONify result.
        result = json.dumps({"rows": result})

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # the cell set to compute differential expression on
        cell_set_base = self.task_def["cellSet"]

        # the cell set we want to compare with (or `rest`)
        cell_set_compare_with = self.task_def["compareWith"]

        # get the top x number of genes to load:
        n_genes = self.task_def.get("maxNum", None)

        # get cell sets from database
        resp = get_item_from_dynamo(self.experiment_id, "cellSets")

        # try to find the right cells
        de_base = find_cells_by_set_id(cell_set_base, resp)

        # use raw values for this task
        raw_adata = self.adata.raw.to_adata()

        if cell_set_compare_with == "rest":
            # We have a simple condition. Everything in the base cluster is `first`,
            # the rest is `second`.
            raw_adata.obs["condition"] = "second"
        else:
            # We have a bit more complicated condition.
            # Everything in the base cluster is `first`,
            # in the second cluster `second`, the rest is `rest`.
            de_compare_with = find_cells_by_set_id(cell_set_compare_with, resp)
            raw_adata.obs["condition"] = "rest"

            raw_adata.obs["condition"].loc[
                raw_adata.obs["cell_ids"].isin(de_compare_with)
            ] = "second"

        raw_adata.obs["condition"].loc[
            raw_adata.obs["cell_ids"].isin(de_base)
        ] = "first"

        # Do a pairwise wilcoxon test
        result = de.test.pairwise(
            data=raw_adata,
            grouping="condition",
            test="rank",
            lazy=False,
            noise_model=None,
        )

        # massage into right format
        result = result.summary_pairs(["first"], ["second"])
        result = result[["gene", "pval", "qval", "log2fc"]]
        result["gene_names"] = result["gene"]
        del result["gene"]

        # remove all NaNs
        result = result.dropna()

        # get top x most significant results, if parameter was supplied
        if n_genes:
            result = result.nsmallest(n_genes, ["qval"])

        return self._format_result(result)
