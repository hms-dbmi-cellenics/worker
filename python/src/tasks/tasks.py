import traceback
import json
from .double_score import DoubletScore
from .mitochondrial_content import MitochondrialContent
from .embedding import ComputeEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
from .cluster_cells import ClusterCells
from result import Result

from config import get_config
from helpers.count_matrix import CountMatrix

config = get_config()


class TaskFactory:
    def __init__(self):
        self.count_matrix = CountMatrix()
        self.count_matrix.sync()

    def submit(self, msg):
        my_class = self._factory(msg)

        # Try to perform task. If fails, send back an error to the API.
        try:
            result = my_class.compute()
            return result
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
            return result

    def _factory(self, msg):
        self.count_matrix.sync()
        adata = self.count_matrix.adata

        if not adata:
            raise Exception("Adata is missing, no tasks can be performed.")

        task_def = msg.get("body", {})
        task_name = task_def.get("name")

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
        elif task_name == "GetDoubletScore":
            my_class = DoubletScore(msg, adata)
            return my_class
        elif task_name == "GetMitochondrialContent":
            my_class = MitochondrialContent(msg, adata)
            return my_class
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
