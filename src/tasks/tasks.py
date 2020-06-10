from .embedding import ComputeEmbedding
from .list_genes import ListGenes
from .differential_expression import DifferentialExpression
from .gene_expression import GeneExpression
import datetime


class TaskFactory:
    def submit(self, msg, adata):

        print(datetime.datetime.now(), "******** ", msg)
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
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
