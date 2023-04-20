import json
import pytest
import responses
from exceptions import PythonWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from worker.tasks.batch_differential_expression import BatchDifferentialExpression
from worker.helpers.s3 import get_cell_sets
from unittest.mock import patch

class TestBatchDifferentialExpression:
    def get_request(
        self,
        cellSet=["cluster1"],
        compareWith="rest",
        basis=["all"],
        comparisonType=None,
    ):
        request = {
            "experimentId": config.EXPERIMENT_ID,
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "BatchDifferentialExpression",
                "cellSet": cellSet,
                "compareWith": compareWith,
                "basis": basis,
            },
        }

        if comparisonType:
            request["body"]["comparisonType"] = comparisonType

        return request

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            BatchDifferentialExpression()

    def test_should_handle_r_worker_error(self):
        basis = ["louvain-0", "louvain-1", "louvain-2", "louvain-3"]
        request_data = self.get_request(
            cellSet=["cluster1"],
            compareWith="cluster2",
            basis=basis,
            comparisonType="within",
        )

        with patch("worker.helpers.s3.get_cell_sets", return_value=cell_set_types["three_sets"]):
            for _ in range(len(basis)):  
                task = BatchDifferentialExpression(request_data)
                result = task.compute()
                assert {"total": 0, "data": "No cell id fulfills the 1st cell set."} in result.data


    def test_should_call_r_worker_endpoint(self):
        cell_sets = ["louvain-0","louvain-1","louvain-2","louvain-3","louvain-4"]
        request_data = self.get_request(
            cellSet=cell_sets,
            compareWith="background",
            basis=["all"],
            comparisonType="within",
        )

        with patch("worker.helpers.s3.get_cell_sets", return_value=cell_set_types["three_sets"]), \
            responses.RequestsMock() as rsps:
                # Add mocked requests for each iteration in the loop
                for _ in range(len(cell_sets)):  # Assuming 4 requests will be made based on the 'basis' value
                    rsps.add(
                        responses.POST,
                        f"{config.R_WORKER_URL}/v0/DifferentialExpression",
                        json={"data": {"full_count": 10, "gene_results": ['ok']}},
                        status=200,
                    )

                task = BatchDifferentialExpression(request_data)
                result = task.compute()
                assert {"total": 10, "data": ['ok']} in result.data