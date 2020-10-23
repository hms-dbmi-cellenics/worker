import pytest
import anndata
import os
from tasks.tasks import TaskFactory
from result import Result
from config import get_config
from mock import Mock, patch

config = get_config()


class TestTaskFactory:
    @pytest.fixture(autouse=True)
    def set_mock_count_matrix_instance(self):
        adata = anndata.read_h5ad(os.path.join(config.LOCAL_DIR, "test", "python.h5ad"))

        with patch("tasks.tasks.CountMatrix") as MockCountMatrix:
            instance = MockCountMatrix.return_value
            instance.sync.return_value = Mock()
            instance.adata.return_value = adata
            self.task_factory = TaskFactory()
            self.task_factory.count_matrix = instance

    def test_throws_typeerror_on_empty_taskfactory_submission(self):
        with pytest.raises(TypeError):
            self.task_factory.submit()

    def test_throws_exception_on_missing_anndata(self):
        self.task_factory.count_matrix.adata = None
        with pytest.raises(Exception) as e:
            self.task_factory.submit({})
        assert e.value.args[0] == "Adata is missing, no tasks can be performed."

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
        assert self.task_factory.count_matrix.sync.call_count == 2

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
        assert self.task_factory.count_matrix.sync.call_count == 2

    def test_each_element_in_result_list_is_a_result_object(self):
        results = self.task_factory.submit(
            {"body": {"name": "GetEmbedding", "type": "pca"}}
        )
        assert all(isinstance(result, Result) for result in results)
        assert self.task_factory.count_matrix.sync.call_count == 2
