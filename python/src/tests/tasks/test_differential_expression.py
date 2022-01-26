import json

import mock
import pytest
import responses
from worker.config import config
from worker.tasks.differential_expression import DifferentialExpression
from worker.helpers.mock_s3 import MockS3Class

class TestDifferentialExpression:
    def get_request(
        self, cellSet="cluster1", compareWith="rest", basis="all", maxNum=None
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

        if maxNum:
            request["body"]["maxNum"] = maxNum

        return request

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        responses.add(
            responses.POST,
            f"{config.R_WORKER_URL}/v0/DifferentialExpression",
            json={},
            status=200,
        )

    """
    Mocks the S3 query for fetching cell sets. Returns an
    empty cell set and yields the patched up object.
    """

    @pytest.fixture
    def mock_S3_get(self):
        with mock.patch("boto3.client") as m:
            mockS3 = MockS3Class()
            m.return_value = mockS3
            yield (m, mockS3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    @responses.activate
    def test_throws_when_second_cellset_missing(self, mock_S3_get):
        MockS3Class.setResponse("one_set")
        with pytest.raises(Exception, match="No cells id fullfills the 2nd cell set"):
            DifferentialExpression(self.get_request())._format_request()

    @responses.activate
    def test_cells_in_sets_intersection_are_filtered_out(self, mock_S3_get):
        MockS3Class.setResponse("two_sets_intersected")

        request = DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="cluster2")
        )._format_request()

        baseCells = request["baseCells"]
        backgroundCells = request["backgroundCells"]


        # Check 1 cell of each of the cell sets is left out
        assert len(baseCells) == len(backgroundCells) == 2

        # Check the cells that haven't been left out are
        # those that are not in the intersection of both sets
        assert len(set(baseCells).intersection(set(backgroundCells))) == 0


    @responses.activate
    def test_cells_not_in_basis_sample_are_filtered_out(self, mock_S3_get):
        MockS3Class.setResponse("three_sets")

        request = DifferentialExpression(
            self.get_request(
                cellSet="cluster1",
                compareWith="cluster2",
                basis="basisCluster",
            )
        )._format_request()

        baseCells = request["baseCells"]
        backgroundCells = request["backgroundCells"]

        # Check cells not in basis are taken out
        assert len(baseCells) == 1
        assert len(backgroundCells) == 2

    @responses.activate
    def test_rest_keyword_only_adds_cells_in_the_same_hierarchy(
        self, mock_S3_get
    ):
        MockS3Class.setResponse("hierarchichal_sets")

        request = DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="rest")
        )._format_request()

        baseCells = request["baseCells"]
        backgroundCells = request["backgroundCells"]

        # Check there is only one cell in each set
        assert len(baseCells) == 1
        assert len(backgroundCells) == 2