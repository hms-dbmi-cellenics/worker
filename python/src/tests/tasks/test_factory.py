import pytest
from mock import Mock, patch

from worker.tasks.factory import TaskFactory


class TestTaskFactory:
    @pytest.fixture(autouse=True)
    def set_mock_task_factory(self):
        with patch("worker.tasks.factory.CountMatrix") as MockCountMatrix:
            instance = MockCountMatrix.return_value
            instance.sync.return_value = Mock()
            self.task_factory = TaskFactory()
            self.task_factory.count_matrix = instance

    def test_throws_typeerror_on_empty_taskfactory_submission(self):
        with pytest.raises(TypeError):
            self.task_factory.submit()

    def test_throws_exception_on_empty_task_definition(self):
        with pytest.raises(Exception) as e:
            self.task_factory.submit({})
        assert e.value.args[0] == "Task class with name None was not found"

    def test_throws_exception_on_non_existent_task(self):
        with pytest.raises(Exception) as e:
            self.task_factory.submit({"body": {"name": "ClearlyAnInvalidTaskName"}})
        assert (
            e.value.args[0]
            == "Task class with name ClearlyAnInvalidTaskName was not found"
        )

    def test_creates_class_on_existent_task(self):
        r = self.task_factory._factory({"body": {"name": "GetEmbedding"}})
        assert isinstance(r, object)
