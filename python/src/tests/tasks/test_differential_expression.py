import json
import os
from unittest.mock import patch

import pytest
import responses
import worker.helpers.s3 as s3
from exceptions import RWorkerException
from worker.config import config
from worker.helpers.mock_s3 import MockS3Class
from worker.tasks.differential_expression import DifferentialExpression


class TestDifferentialExpression:
    def get_request(
        self,
        cellSet="cluster1",
        compareWith="rest",
        basis="all",
        comparisonType=None,
        maxNum=None,
    ):
        request = {
            "experimentId": "e52b39624588791a7889e39c617f669e1",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "DifferentialExpression",
                "cellSet": cellSet,
                "compareWith": compareWith,
                "basis": basis,
            },
        }

        if comparisonType:
            request["body"]["comparisonType"] = comparisonType

        if maxNum:
            request["body"]["maxNum"] = maxNum

        return request

    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "MockCellSet.json")) as f:
            self.cellsets = json.load(f)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    # @responses.activate
    # def test_cells_not_in_basis_sample_are_filtered_out(self, mock_S3_get):
    #     MockS3Class.setResponse("three_sets")

    #     request = DifferentialExpression(
    #         self.get_request(
    #             cellSet="cluster1",
    #             compareWith="cluster2",
    #             basis="basisCluster",
    #         )
    #     )._format_request()

    #     baseCells = request["baseCells"]
    #     backgroundCells = request["backgroundCells"]

    #     # Check cells not in basis are taken out
    #     assert len(baseCells) == 1
    #     assert len(backgroundCells) == 2

    # @responses.activate
    # def test_rest_keyword_only_adds_cells_in_the_same_hierarchy(self, mock_S3_get):
    #     MockS3Class.setResponse("hierarchichal_sets")

    #     request = DifferentialExpression(
    #         self.get_request(cellSet="cluster1", compareWith="rest")
    #     )._format_request()

    #     baseCells = request["baseCells"]
    #     backgroundCells = request["backgroundCells"]

    #     # Check there is only one cell in each set
    #     assert len(baseCells) == 1
    #     assert len(backgroundCells) == 2

    # @responses.activate
    # def test_default_comparison_type_added_to_request(self, mock_S3_get):

    #     request = DifferentialExpression(
    #         self.get_request(
    #             cellSet="cluster1",
    #             compareWith="cluster2",
    #             basis="basisCluster",
    #         )
    #     )._format_request()

    #     # Check that comparisonType defaults to within
    #     comparisonType = request["comparisonType"]
    #     assert comparisonType == "within"

    @responses.activate
    @patch("worker.helpers.s3.get_cell_sets")
    def test_should_throw_exception_on_r_worker_error(self, mock_get_cell_sets):
        error_code = "R_WORKER_ERROR"
        user_message = "User message"

        mock_get_cell_sets.return_value = self.cellsets

        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            json={
                "error": {
                    "error_code": error_code,
                    "user_message": user_message,
                }
            },
            status=200,
        )

        with pytest.raises(RWorkerException) as exc_info:
            DifferentialExpression(
                self.get_request(
                    cellSet="louvain-1",
                    compareWith="louvain-2",
                    basis="condition-control",
                )
            ).compute()

        assert exc_info.value.args[0] == error_code
        assert exc_info.value.args[1] == user_message
