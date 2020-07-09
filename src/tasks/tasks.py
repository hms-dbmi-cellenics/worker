from .embedding import ComputeEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
from .prepare_experiment import PrepareExperiment

from helpers import load_count_matrix


class TaskFactory:
    def submit(self, msg, adata):
        try:
            my_class = self._factory(msg, adata)
            result = my_class.compute()
            return result
        except Exception as e:
            # do return this though to the api
            raise e

    @staticmethod
    def _factory(msg, adata):
        task_def = msg["body"]
        task_name = task_def["name"]

        if not adata and task_name != "PrepareExperiment":
            adata = load_count_matrix.get_adata(adata, msg["experimentId"])

        if task_name == "GetEmbedding":
            my_class = ComputeEmbedding(msg, adata)
            return my_class
        if task_name == "ListGenes":
            my_class = ListGenes(msg, adata)
            return my_class
        if task_name == "DifferentialExpression":
            my_class = DifferentialExpression(msg, adata)
            return my_class
        if task_name == "GeneExpression":
            my_class = GeneExpression(msg, adata)
            return my_class
        if task_name == "PrepareExperiment":
            # after this line, adata will be equal to the newly computed adata file
            my_class = PrepareExperiment(msg, adata)
            return my_class
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
