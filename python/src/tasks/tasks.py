import traceback
import json
import anndata
from .embedding import ComputeEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
from .cluster_cells import ClusterCells
from result import Result

from config import get_config
from helpers import count_matrix

config = get_config()


class TaskFactory:
    def __init__(self, experimentId):
        self.adata = None
        self.adata_path = count_matrix.get_adata_path()
        self._initialise_adata()

    def _initialise_adata(self):
        count_matrix.download_all_files()
        with open(self.adata_path, "rb+") as f:
            self.adata = anndata.read_h5ad(f)
            if "cell_ids" not in self.adata.obs:
                raise ValueError(
                    "You must have `cell_ids` in your anndata file for integer cell IDs."
                )

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
        print("before checking for change: ", self.adata_path)
        if count_matrix.is_file_changed(self.adata, self.adata_path):
            print("There has been recent write to the file, have to download it again.")
            # TODO: after introducing multiple-sample support, make this
            #  more efficient by downloading only a specific file that got changed.
            self._initialise_adata()
        else:
            print("The Anndata file is the same as before, no need to download it.")

        task_def = msg["body"]
        task_name = task_def["name"]

        if task_name == "GetEmbedding":
            my_class = ComputeEmbedding(msg, self.adata)
            return my_class
        elif task_name == "ListGenes":
            my_class = ListGenes(msg, self.adata)
            return my_class
        elif task_name == "DifferentialExpression":
            my_class = DifferentialExpression(msg, self.adata)
            return my_class
        elif task_name == "GeneExpression":
            my_class = GeneExpression(msg, self.adata)
            return my_class
        elif task_name == "ClusterCells":
            my_class = ClusterCells(msg, self.adata)
            return my_class
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
