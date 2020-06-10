import pandas as pd
import numpy as np
import json
import scanpy
import boto3
from config import get_config
from result import Result

config = get_config()


class DifferentialExpression:
    def __init__(self, msg, adata):
        self.adata = adata
        self.dynamo = boto3.resource("dynamodb").Table(config.get_dynamo_table())

        self.task_def = msg["body"]
        self.experiment_id = msg["experimentId"]

    def find_cells_by_set_id(self, needle, haystack):
        for cell_set in haystack:
            if cell_set["key"] == needle:
                return cell_set["cellIds"]

            if cell_set.get("children", None):
                result = self.find_cells_by_set_id(needle, cell_set["children"])

                if result:
                    return result

        return []

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
        MAX_NUM_GENES = 10e11
        n_genes = self.task_def.get("maxNum", MAX_NUM_GENES)

        # get cell sets from database
        resp = self.dynamo.get_item(
            Key={"experimentId": self.experiment_id}, ProjectionExpression="cellSets",
        )
        resp = resp["Item"]["cellSets"]

        # try to find the right cells
        de_base = self.find_cells_by_set_id(cell_set_base, resp)
        self.adata.obs["de_base"] = self.adata.obs.index.isin(de_base)
        self.adata.obs["de_base"] = pd.Categorical(
            self.adata.obs["de_base"], categories=[True]
        )

        # if `compareWith` is not 'rest', try create another group based on that set
        if cell_set_compare_with != "rest":
            de_compare_with = self.find_cells_by_set_id(cell_set_compare_with, resp)
            self.adata.obs["de_compare_with"] = self.adata.obs.index.isin(
                de_compare_with
            )

            self.adata.obs["de_compare_with"] = pd.Categorical(
                self.adata.obs["de_compare_with"], categories=[True]
            )

        scanpy.tl.rank_genes_groups(
            self.adata, "de_base", method="t-test", n_genes=n_genes
        )
        de_result = self.adata.uns["rank_genes_groups"]

        de_result["gene_names"] = de_result["names"]
        del de_result["names"]
        del de_result["params"]

        de_result = pd.DataFrame.from_dict(de_result)
        de_result = de_result.apply(pd.Series.explode)

        print(de_result)

        return self._format_result(de_result)
