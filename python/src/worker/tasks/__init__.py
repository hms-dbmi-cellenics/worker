from abc import ABC, abstractmethod


class Task(ABC):
    """A task submitted to the worker."""

    def __init__(self, msg):
        self.task_def = msg["body"]

    @abstractmethod
    def compute(self):
        ...

    @abstractmethod
    def _format_result(self, result):
        ...


RDS_PATH = "/data/RResults/r.rds"