import pytest
import os
from tasks.tasks import TaskFactory
from result import Result
from config import get_config
from mock import Mock, patch
from helpers.count_matrix import CountMatrix

config = get_config()


class TestTaskFactory:
    @pytest.fixture(autouse=True)
    def set_mock_task_factory(self):
        with patch("tasks.tasks.CountMatrix") as MockCountMatrix:
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

    @pytest.mark.parametrize(
        "valid_task_name",
        [
            "GetEmbedding",
            "ListGenes",
            "DifferentialExpression",
            "GeneExpression",
            "ClusterCells",
        ],
    )
    def test_returns_result_list_with_properly_defined_task(self, valid_task_name):
        results = self.task_factory.submit({"body": {"name": valid_task_name}})
        assert isinstance(results, list)

    def test_each_element_in_result_list_is_a_result_object(self):
        results = self.task_factory.submit(
            {"body": {"name": "GetEmbedding", "type": "pca"}}
        )
        assert all(isinstance(result, Result) for result in results)
