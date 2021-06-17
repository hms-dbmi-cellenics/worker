from abc import ABC, abstractmethod


class Task(ABC):
    """A task submitted to the worker."""
    def __init__(self, msg=None):
        pass

    @abstractmethod
    def compute(self):
        ...

    @abstractmethod
    def _format_result(self, result):
        ...
