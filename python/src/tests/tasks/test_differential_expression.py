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

    def get_s3_stub(self, content_type):
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
        with open(os.path.join("tests/data", "cell_set_types.json")) as f:
            all_cell_set_types = json.load(f)
        content = all_cell_set_types[content_type]
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
        return (stubber, s3)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    def test_throws_when_second_cellset_missing(self):
        stubber, s3 = self.get_s3_stub("one_set")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            with pytest.raises(
                Exception, match="No cell id fullfills the 2nd cell set"
            ):
                DifferentialExpression(self.get_request())._format_request()
            stubber.assert_no_pending_responses()

    def test_cells_in_sets_intersection_are_filtered_out(self):
        stubber, s3 = self.get_s3_stub("two_sets_intersected")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
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
            stubber.assert_no_pending_responses()

    def test_cells_not_in_basis_sample_are_filtered_out(self):
        stubber, s3 = self.get_s3_stub("three_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
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
            stubber.assert_no_pending_responses()

    def test_rest_keyword_only_adds_cells_in_the_same_hierarchy(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = DifferentialExpression(
                self.get_request(cellSet="cluster1", compareWith="rest")
            )._format_request()

            baseCells = request["baseCells"]
            backgroundCells = request["backgroundCells"]

            # Check there is only one cell in each set
            assert len(baseCells) == 1
            assert len(backgroundCells) == 2
            stubber.assert_no_pending_responses()

    def test_default_comparison_type_added_to_request(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = DifferentialExpression(
                self.get_request(
                    cellSet="cluster1",
                    compareWith="cluster2",
                    basis="basisCluster",
                )
            )._format_request()

            # Check that comparisonType defaults to within
            comparisonType = request["comparisonType"]
            assert comparisonType == "within"
            stubber.assert_no_pending_responses()

    def test_specified_comparison_type_added_to_request(self):
        stubber, s3 = self.get_s3_stub("hierarchichal_sets")

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = DifferentialExpression(
                self.get_request(
                    cellSet="cluster1",
                    compareWith="cluster2",
                    basis="basisCluster",
                    comparisonType="between",
                )
            )._format_request()

            # Check that comparisonType uses set value of between instead of default (within)
            comparisonType = request["comparisonType"]
            assert comparisonType == "between"
            stubber.assert_no_pending_responses()
