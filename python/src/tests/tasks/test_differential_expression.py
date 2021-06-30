import json
from operator import itemgetter

import mock
import pytest
import responses

from worker.config import config
from worker.tasks.differential_expression import DifferentialExpression

cell_set_responses = {
    "one_set": [
        {"name": "my amazing cluster", "key": "cluster1", "cellIds": [4, 5]},
        {
            "name": "my other amazing cluster",
            "key": "cluster2",
            "cellIds": [0, 1, 2, 3],
        },
    ],
    "two_sets": [
        {
            "name": "my amazing cluster",
            "key": "cluster1",
            "cellIds": [4, 5],
        },
        {
            "name": "my other amazing cluster",
            "key": "cluster2",
            "cellIds": [0, 1, 2, 3],
        },
    ],
    "two_sets_intersected": [
        {
            "name": "intersecting set",
            "key": "cluster1",
            "cellIds": [1, 2, 3],
        },
        {
            "name": "other intersecting set",
            "key": "cluster2",
            "cellIds": [3, 4, 5],
        },
    ],
    "three_sets": [
        {
            "name": "one set",
            "key": "cluster1",
            "cellIds": [4, 5],
        },
        {
            "name": "other set",
            "key": "cluster2",
            "cellIds": [0, 1, 2, 3],
        },
        {
            "name": "basis set",
            "key": "basisCluster",
            "cellIds": [0, 1, 5],
        },
    ],
    "hierarchichal_sets": [
        {
            "name": "hierarchy 1",
            "key": "set_hierarchy_1",
            "cellIds": [],
            "children": [
                {
                    "name": "one set",
                    "key": "cluster1",
                    "cellIds": [4],
                },
                {
                    "name": "another set",
                    "key": "cluster2",
                    "cellIds": [5],
                },
            ],
        },
        {
            "name": "hierarchy 2",
            "key": "set_hierarchy_2",
            "cellIds": [],
            "children": [
                {
                    "name": "set",
                    "key": "cluster3",
                    "cellIds": [0],
                },
                {
                    "name": "set1",
                    "key": "cluster4",
                    "cellIds": [1],
                },
            ],
        },
    ],
}


class MockS3Class:
    response = cell_set_responses["one_set"]

    no_called = 0

    def setResponse(response_key):
        MockS3Class.response = cell_set_responses[response_key]

    # def Table(*args, **kwargs):
    #     MockS3Class.no_called = 0
    #     return MockS3Class.MockTable()

    # def download_fileobj(*args, **kwargs):
    #     MockS3Class.no_called += 1
    #     return {"Body": json.dumps({"cellSets": MockS3Class.response})}

    def download_fileobj(*args, **kwargs):
        Bucket, Key, Fileobj = itemgetter("Bucket", "Key", "Fileobj")(kwargs)

        if not Bucket or not Key or not Fileobj:
            raise Exception("Parameters not received")

        Fileobj.write(str.encode(f'{{"cellSets": {json.dumps(MockS3Class.response)}}}'))

        return
        # MockS3Class.no_called += 1
        # return {"Body": json.dumps({"cellSets": MockS3Class.response})}


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
    def test_cell_sets_get_queried_appropriately(self, mock_S3_get):
        DifferentialExpression(self.get_request()).compute()

    @responses.activate
    def test_works_when_all_is_first(self, mock_S3_get):
        request = self.get_request(cellSet="all-asdasd", compareWith="cluster1")

        DifferentialExpression(request)

    @responses.activate
    def test_returns_json(self, mock_S3_get):
        res = DifferentialExpression(self.get_request()).compute()
        res = res[0].result
        json.loads(res)

    @responses.activate
    def test_returns_a_json_object(self, mock_S3_get):
        res = DifferentialExpression(self.get_request()).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    @responses.activate
    def test_object_has_all_required_columns(self, mock_S3_get):
        res = DifferentialExpression(self.get_request()).compute()
        res = res[0].result
        res = json.loads(res)

        for row in res["rows"]:
            keys = sorted(row.keys())
            expected_keys = sorted(
                # Until the UI side is not changed we need to suppor old and new columns
                # ["gene_names", "_row", "avg_log2FC", "p_val_adj", "pct_1", "pct_2"]
                [
                    "gene_names",
                    "zscore",
                    "abszscore",
                    "qval",
                    "log2fc",
                    "_row",
                    "avg_log2FC",
                    "p_val_adj",
                    "pct_1",
                    "pct_2",
                ]
            )
            print(expected_keys)
            print(keys)
            assert keys == expected_keys

    @responses.activate
    def test_cells_in_sets_intersection_are_filtered_out(self, mock_S3_get):
        MockS3Class.setResponse("two_sets_intersected")

        DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="cluster2")
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        print(baseCells)
        print(backgroundCells)

        # Check 1 cell of each of the cell sets is left out
        assert len(baseCells) == len(backgroundCells) == 2

        # Check the cells that haven't been left out are
        # those that are not in the intersection of both sets
        assert len(set(baseCells).intersection(set(backgroundCells))) == 0

    @responses.activate
    def test_cells_not_in_basis_sample_are_filtered_out(self, mock_S3_get):
        MockS3Class.setResponse("three_sets")

        DifferentialExpression(
            self.get_request(
                cellSet="cluster1", compareWith="cluster2", basis="basisCluster"
            )
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        # Check cells not in basis are taken out
        assert len(baseCells) == 1
        assert len(backgroundCells) == 2

    @responses.activate
    def test_rest_keyword_only_adds_cells_in_the_same_hierarchy(self, mock_S3_get):
        MockS3Class.setResponse("hierarchichal_sets")

        DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="rest")
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        # Check there is only one cell in each set
        assert len(baseCells) == 1
        assert len(backgroundCells) == 1


"""
    def test_dynamodb_call_is_made_once_when_vs_rest(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        MockDynamoClass.setResponse("two_sets")

        DifferentialExpression(self.get_request()).compute()

        assert dynamodb.no_called == 1

    def test_appropriate_genes_returned_when_a_limit_is_specified(
        self, mock_dynamo_get
    ):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        request = self.get_request(maxNum=2)

        res = DifferentialExpression(request).compute()
        res = res[0].result
        res = json.loads(res)["rows"]

        assert len(res) <= request["body"]["maxNum"]

    # In these three tests we don't actually care about the end result of the r worker, we just need to see the request generated
    # on the python side, so we can leave responses.activate enabled with an empty response {}.
"""