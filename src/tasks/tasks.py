from .embedding import ComputeEmbedding


class TaskFactory:
    def factory(task_type, *args, **kwargs):
        # return eval(type + "()")
        if task_type == "ComputeEmbedding":
            return ComputeEmbedding(*args, **kwargs)
        else:
            raise Exception("Task class with name {} was not found".format(task_type))
