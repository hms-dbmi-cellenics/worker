import io
import json
import os

import boto3
import mock
import pytest
from botocore.stub import Stubber

from worker.config import config
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
            "experimentId": config.EXPERIMENT_ID,
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

    """
    Mocks the S3 query for fetching cell sets. Returns an
    empty cell set and yields the patched up object.
    """

    # @pytest.fixture
    # def mock_S3_get(self):

    #     return stubber

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    def test_throws_when_second_cellset_missing(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        response = {
            "ContentLength": 10,
            "ContentType": "utf-8",
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }

        expected_params = {
            "Bucket": config.CELL_SETS_BUCKET,
            "Key": config.EXPERIMENT_ID,
        }
        stubber = Stubber(s3)
        stubber.add_response("head_object", response, expected_params)

        # Get object
        content = {
            "cellSets": [
                {"name": "my amazing cluster", "key": "cluster1", "cellIds": [4, 5]}
            ]
        }

        content_bytes = json.dumps(content, indent=2).encode("utf-8")
        # with open(os.path.join("tests/data", "cell_sets.json"), "rb") as f:
        #     content = f.read()
        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        response = {
            "ContentLength": len(content_bytes),
            "ContentType": "utf-8",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.CELL_SETS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            with pytest.raises(
                Exception, match="No cell id fullfills the 2nd cell set"
            ):
                DifferentialExpression(self.get_request())._format_request()

    # def test_cells_in_sets_intersection_are_filtered_out(self, mock_S3_get):
    #     MockS3Class.setResponse("two_sets_intersected")

    #     request = DifferentialExpression(
    #         self.get_request(cellSet="cluster1", compareWith="cluster2")
    #     )._format_request()

    #     baseCells = request["baseCells"]
    #     backgroundCells = request["backgroundCells"]

    #     # Check 1 cell of each of the cell sets is left out
    #     assert len(baseCells) == len(backgroundCells) == 2

    #     # Check the cells that haven't been left out are
    #     # those that are not in the intersection of both sets
    #     assert len(set(baseCells).intersection(set(backgroundCells))) == 0

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

    # def test_specified_comparison_type_added_to_request(self, mock_S3_get):

    #     request = DifferentialExpression(
    #         self.get_request(
    #             cellSet="cluster1",
    #             compareWith="cluster2",
    #             basis="basisCluster",
    #             comparisonType="between",
    #         )
    #     )._format_request()

    #     # Check that comparisonType uses set value of between instead of default (within)
    #     comparisonType = request["comparisonType"]
    #     assert comparisonType == "between"
