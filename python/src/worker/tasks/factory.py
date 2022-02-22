import json
import traceback
from logging import error

from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.count_matrix import CountMatrix
from ..helpers.worker_exception import WorkerException
from ..result import Result
from ..tasks import Task
from .background_expressed_genes import GetBackgroundExpressedGenes
from .cluster_cells import ClusterCells
from .differential_expression import DifferentialExpression
from .dotplot import DotPlot
from .doublet_score import GetDoubletScore
from .embedding import GetEmbedding
from .expression_cellsets import GetExpressionCellSets
from .gene_expression import GeneExpression
from .list_genes import ListGenes
from .marker_heatmap import MarkerHeatmap
from .mitochondrial_content import GetMitochondrialContent


class TaskFactory:
    tasks = {
        t.__name__: t
        for t in (
            GetEmbedding,
            ListGenes,
            DifferentialExpression,
            GeneExpression,
            GetBackgroundExpressedGenes,
            ClusterCells,
            DotPlot,
            GetDoubletScore,
            GetMitochondrialContent,
            MarkerHeatmap,
            GetExpressionCellSets,
        )
    }

    def __init__(self):
        self.count_matrix = CountMatrix()
        self.count_matrix.sync()

    def _log_exception(self, task, error_object):
        trace = (
            f"Exception for task {task.__class__.__name__}:\n"
            f"{traceback.format_exc()}"
        )

        error(trace)
        xray_recorder.current_segment().add_exception(error_object, trace)

        return trace

    def submit(self, msg):
        task = self._factory(msg)

        # Try to perform task. If fails, send back an error to the API.
        try:
            result = task.compute()
            return result

        except Exception as e:
            trace = self._log_exception(task, e)

            if config.CLUSTER_ENV == "development":
                result = Result(trace, error=True)

            if isinstance(e, WorkerException):
                return Result(
                    {
                        "error_code": e.error_code,
                        "user_message": e.user_message,
                    },
                    error=True,
                )

                # Only send real traces in development.

            return Result(
                {
                    "error_code": "PYTHON_WORKER_GENERIC_ERROR",
                    "user_message": "An unexpected error occurred while performing the work.",
                },
                error=True,
            )

    def _factory(self, msg) -> Task:
        self.count_matrix.sync()
        task_def = msg.get("body", {})
        task_name = task_def.get("name")

        try:
            return self.tasks[task_name](msg)
        except KeyError as e:
            raise ValueError(f"Task class with name {task_name} was not found") from e
