import pytest
import anndata
import os
from tasks.tasks import TaskFactory
from result import Result


class TestResult:
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

    def test_throws_exception_on_incomplete_name(self):
        with pytest.raises(KeyError):
            TaskFactory().submit({"type": "asd"}, self._adata)

    def test_throws_exception_on_incomplete_type(self):
        with pytest.raises(KeyError):
            TaskFactory().submit({"name": "asd"}, self._adata)

    def test_throws_exception_on_non_existent_task(self):
        with pytest.raises(Exception):
            TaskFactory().submit(
                {"name": "ClearlyAnInvalidTaskName", "type": "asd"}, self._adata
            )

    def test_creates_class_on_existent_task(self):
        r = TaskFactory._factory("GetEmbedding", self._adata)
        assert isinstance(r, object)

    def test_returns_result_list_with_properly_defined_task(self):
        results = TaskFactory().submit(
            {"name": "GetEmbedding", "type": "pca"}, self._adata
        )
        assert isinstance(results, list)

    def test_each_element_in_result_list_is_a_result_object(self):
        results = TaskFactory().submit(
            {"name": "GetEmbedding", "type": "pca"}, self._adata
        )
        assert all(isinstance(result, Result) for result in results)
