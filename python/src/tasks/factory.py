import traceback
import json
from logging import error

from .doublet_score import GetDoubletScore
from .mitochondrial_content import GetMitochondrialContent
from .embedding import GetEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
from .cluster_cells import ClusterCells
from result import Result
from tasks import Task
from aws_xray_sdk.core import xray_recorder

from config import config
from helpers.count_matrix import CountMatrix


class TaskFactory:
    tasks = {t.__name__: t for t in (GetEmbedding,
                                     ListGenes,
                                     DifferentialExpression,
                                     GeneExpression,
                                     ClusterCells,
                                     GetDoubletScore,
                                     GetMitochondrialContent)
             }

    def __init__(self):
        self.count_matrix = CountMatrix()
        self.count_matrix.sync()

    def submit(self, msg):
        task = self._factory(msg)

        # Try to perform task. If fails, send back an error to the API.
        try:
            result = task.compute()
            return result
        except Exception as e:
            trace = traceback.format_exc()
            error(trace)

            xray_recorder.current_segment().add_exception(e, trace)

            # Only send real traces in development.
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

    def _factory(self, msg) -> Task:
        self.count_matrix.sync()
        task_def = msg.get("body", {})
        task_name = task_def.get("name")

        try:
            return self.tasks[task_name]
        except KeyError as e:
            raise ValueError(f"Task class with name {task_name} was not found") from e
