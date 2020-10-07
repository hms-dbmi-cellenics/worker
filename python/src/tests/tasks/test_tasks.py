import pytest
import anndata
import os
from tasks.tasks import TaskFactory
from result import Result

from helpers import load_count_matrix
from mock import Mock


class TestTaskFactory:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(os.path.join("tests", "test.h5ad"))

    def test_throws_typeerror_on_empty_taskfactory_submission(self):
        with pytest.raises(TypeError):
            TaskFactory().submit()

    def test_throws_exception_on_missing_anndata(self):
        with pytest.raises(TypeError):
            TaskFactory().submit({})

    def test_throws_exception_on_empty_task_definition(self):
        with pytest.raises(TypeError):
            TaskFactory().submit({})

    def test_throws_exception_on_incomplete_body(self):
        with pytest.raises(KeyError):
            TaskFactory().submit({"body": {}}, self._adata)

    def test_throws_exception_on_non_existent_task(self):
        with pytest.raises(Exception):
            TaskFactory().submit(
                {"body": {"name": "ClearlyAnInvalidTaskName"}}, self._adata
            )

    def test_creates_class_on_existent_task(self):
        r = TaskFactory._factory({"body": {"name": "GetEmbedding"}}, self._adata)
        assert isinstance(r, object)

    def test_returns_result_list_with_properly_defined_task(self):
        results, adata = TaskFactory().submit(
            {"body": {"name": "GetEmbedding", "type": "pca"}}, self._adata
        )
        assert isinstance(results, list)

    def test_each_element_in_result_list_is_a_result_object(self):
        results, adata = TaskFactory().submit(
            {"body": {"name": "GetEmbedding", "type": "pca"}}, self._adata
        )

        assert all(isinstance(result, Result) for result in results)

    def loads_adata_if_not_loaded(self):
        load_count_matrix.get_adata = Mock(return_value=self._adata)

        TaskFactory._factory(
            {"body": {"name": "GetEmbedding"}, "experimentId": "1234"}, None
        )
        load_count_matrix.get_adata.assert_called_once_with(None, "1234")

    def test_dont_load_adata_if_loaded(self):
        load_count_matrix.get_adata = Mock()

        TaskFactory._factory(
            {"body": {"name": "GetEmbedding"}, "experimentId": "1234"}, self._adata
        )
        assert not load_count_matrix.get_adata.called

    def test_dont_load_adata_if_not_needed(self):
        load_count_matrix.get_adata = Mock()

        TaskFactory._factory(
            {
                "body": {
                    "name": "PrepareExperiment",
                    "sourceMatrixPath": "my/path",
                    "sourceBucket": "my-bucket",
                },
                "experimentId": "1234",
            },
            None,
        )
        assert not load_count_matrix.get_adata.called
