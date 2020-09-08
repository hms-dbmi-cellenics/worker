import traceback
import json

from .embedding import ComputeEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
from .prepare_experiment import PrepareExperiment
from .cluster_cells import ClusterCells
from result import Result

from config import get_config
from helpers import load_count_matrix

config = get_config()


class TaskFactory:
    def submit(self, msg, adata):
        if not adata and msg["body"]["name"] != "PrepareExperiment":
            adata = load_count_matrix.get_adata(adata, msg["experimentId"])
        my_class = self._factory(msg, adata)

        # Try to perform task. If fails, send back an error to the API.
        try:
            result = my_class.compute()
            return result, adata
        except Exception:
            trace = traceback.format_exc()
            print(trace)

            # Do not send real traces in development.
            if config.CLUSTER_ENV == "development":
                result = [Result(json.dumps(trace), error=True)]
            else:
                result = [
                    Result(
                        json.dumps(
                            "An unexpected error occurred while performing the work."
                        ),
                        error=True,
                    )
                ]
            return result, adata

    @staticmethod
    def _factory(msg, adata):
        task_def = msg["body"]
        task_name = task_def["name"]

        if task_name == "GetEmbedding":
            my_class = ComputeEmbedding(msg, adata)
            return my_class
        elif task_name == "ListGenes":
            my_class = ListGenes(msg, adata)
            return my_class
        elif task_name == "DifferentialExpression":
            my_class = DifferentialExpression(msg, adata)
            return my_class
        elif task_name == "GeneExpression":
            my_class = GeneExpression(msg, adata)
            return my_class
        elif task_name == "ClusterCells":
            my_class = ClusterCells(msg, adata)
            return my_class
        elif task_name == "PrepareExperiment":
            # after this line, adata will be equal to the newly computed adata file
            my_class = PrepareExperiment(msg, adata)
            return my_class
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
