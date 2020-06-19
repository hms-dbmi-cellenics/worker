import json
import scanpy
import boto3
from config import get_config
from result import Result

from tasks.helpers.find_cells_by_set_id import find_cells_by_set_id

config = get_config()


class GeneExpression:
    def __init__(self, msg, adata):
        self.adata = adata
        self.dynamo = boto3.resource("dynamodb", region_name=config.AWS_REGION).Table(
            config.get_dynamo_table()
        )

        self.task_def = msg["body"]
        self.experiment_id = msg["experimentId"]

    def _aggregate_cells_from_cell_sets(self, cell_sets):
        # get cell sets from database
        resp = self.dynamo.get_item(
            Key={"experimentId": self.experiment_id}, ProjectionExpression="cellSets",
        )
        resp = resp["Item"]["cellSets"]

        cells = set()

        for cell_set in cell_sets:
            cells_found = find_cells_by_set_id(cell_set, resp)
            cells.update(cells_found)

        return cells

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)

        # Return a list of formatted results.
        return [Result(result)]

    def compute(self):
        # the cell sets to get expression data from
        cell_sets = self.task_def.get("cellSets", [])

        # the genes to get expression data for
        genes = self.task_def["genes"]

        # whether to perform feature scaling (defaults to False)
        scale = self.task_def.get("scale", False)

        raw_adata = self.adata.raw.to_adata()

        # try to find all cells in the list
        if cell_sets == "all":
            raw_adata.obs["cells_to_compute"] = True
        else:
            cells = self._aggregate_cells_from_cell_sets(cell_sets)
            raw_adata.obs["cells_to_compute"] = raw_adata.obs.index.isin(cells)

        # try to find all genes in the list
        raw_adata.var["genes_to_compute"] = raw_adata.var.index.isin(genes)

        # filter matrix
        raw_adata = raw_adata[
            raw_adata.obs["cells_to_compute"], raw_adata.var["genes_to_compute"]
        ]

        # if feature scaling is desired, perform that now
        if scale:
            scanpy.pp.scale(raw_adata, max_value=10)

        # compute result
        result = {"cells": raw_adata.obs.index.tolist(), "data": []}

        for gene in genes:
            view = raw_adata[:, raw_adata.var.index == gene]
            expression = view.X.toarray().flatten().tolist()

            result["data"].append({"geneName": gene, "expression": expression})

        return self._format_result(result)
