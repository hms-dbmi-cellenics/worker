from exceptions import WorkerException
from worker_status_codes import PYTHON_WORKER_ERROR

from ..helpers.count_matrix import CountMatrix
from ..helpers.xray_log_exception import xray_log_exception
from ..result import Result
from ..tasks import Task
from .background_expressed_genes import GetBackgroundExpressedGenes
from .cluster_cells import ClusterCells
from .differential_expression import DifferentialExpression
from .batch_differential_expression import BatchDifferentialExpression
from .dotplot import DotPlot
from .doublet_score import GetDoubletScore
from .embedding import GetEmbedding
from .expression_cellsets import GetExpressionCellSets
from .gene_expression import GeneExpression
from .list_genes import ListGenes
from .marker_heatmap import MarkerHeatmap
from .mitochondrial_content import GetMitochondrialContent
from .n_genes import GetNGenes
from .n_umis import GetNUmis
from .trajectory_analysis_pseudotime import GetTrajectoryAnalysisPseudoTime
from .trajectory_analysis_starting_nodes import GetTrajectoryAnalysisStartingNodes
from .normalized_matrix import GetNormalizedExpression
from .cell_annotation_sctype import ScTypeAnnotate
from .download_annot_seurat_object import DownloadAnnotSeuratObject


class TaskFactory:
    tasks = {
        t.__name__: t
        for t in (
            GetEmbedding,
            ListGenes,
            DifferentialExpression,
            BatchDifferentialExpression,
            GeneExpression,
            GetBackgroundExpressedGenes,
            ClusterCells,
            DotPlot,
            GetDoubletScore,
            GetMitochondrialContent,
            GetNGenes,
            GetNUmis,
            MarkerHeatmap,
            GetTrajectoryAnalysisStartingNodes,
            GetTrajectoryAnalysisPseudoTime,
            GetExpressionCellSets,
            GetNormalizedExpression,
            ScTypeAnnotate,
            DownloadAnnotSeuratObject
        )
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

        except WorkerException as e:
            xray_log_exception(task, e)
            return Result(
                {
                    "error_code": e.error_code,
                    "user_message": e.user_message,
                },
                error=True,
            )

        except Exception as e:
            xray_log_exception(task, e)
            return Result(
                {
                    "error_code": PYTHON_WORKER_ERROR,
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
            raise KeyError(f"Task class with name {task_name} was not found: {e}")
