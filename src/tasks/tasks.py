from .embedding import ComputeEmbedding
import datetime


class TaskFactory:
    def submit(self, task_def, adata):
        task_name = task_def["name"]
        details = task_def["type"]

        print(datetime.datetime.now(), "******** ", task_name)
        try:
            my_class = self._factory(task_name, adata)
            result = my_class.compute(details)
            return result
        except Exception as e:
            # do return this though to the api
            raise e

    @staticmethod
    def _factory(task_name, adata):
        # return eval(type + "()")
        if task_name == "GetEmbedding":
            my_class = ComputeEmbedding(adata)
            return my_class
        else:
            raise Exception("Task class with name {} was not found".format(task_name))
