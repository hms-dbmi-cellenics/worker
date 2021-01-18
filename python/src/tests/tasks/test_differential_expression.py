import pytest
import anndata
import os
from tasks.differential_expression import DifferentialExpression
import json
import mock
import responses
from config import get_config

config = get_config()

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


class MockDynamoClass:
    response = cell_set_responses["one_set"]

    no_called = 0

    def setResponse(response_key):
        MockDynamoClass.response = cell_set_responses[response_key]

    def Table(*args, **kwargs):
        MockDynamoClass.no_called = 0
        return MockDynamoClass.MockTable()

    class MockTable:
        def get_item(*args, **kwargs):
            MockDynamoClass.no_called += 1
            print("doing stuff")
            return {"Item": {"cellSets": MockDynamoClass.response}}


class TestDifferentialExpression:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(
            os.path.join(config.LOCAL_DIR, "test", "python.h5ad")
        )

    def get_request(
        self, cellSet="cluster1", compareWith="rest", basis="all", maxNum=None
    ):
        request = {
            "experimentId": "5e959f9c9f4b120771249001",
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
        with open(os.path.join("tests", "de_result.json")) as f:
            data = json.load(f)
            responses.add(
                responses.POST,
                f"{config.R_WORKER_URL}/v0/DifferentialExpression",
                json=data,
                status=200,
            )

    """
    Mocks the DynamoDB query for fetching cell sets. Returns an
    empty cell set and yields the patched up object.
    """

    @pytest.fixture
    def mock_dynamo_get(self):
        with mock.patch("boto3.resource") as m:
            mockDynamo = MockDynamoClass()
            m.return_value = mockDynamo
            yield (m, mockDynamo)

    @responses.activate
    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DifferentialExpression()

    @responses.activate
    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            DifferentialExpression(self.get_request())

    @responses.activate
    def test_dynamodb_call_is_made_once_when_vs_rest(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        MockDynamoClass.setResponse("two_sets")

        DifferentialExpression(self.get_request(), self._adata).compute()

        assert dynamodb.no_called == 1

    @responses.activate
    def test_cell_sets_get_queried_appropriately(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb
        DifferentialExpression(self.get_request(), self._adata).compute()

    @responses.activate
    def test_works_with_request_and_adata(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb
        DifferentialExpression(self.get_request(), self._adata)

    @responses.activate
    def test_works_when_all_is_first(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        request = self.get_request(cellSet="all-asdasd", compareWith="cluster1")

        DifferentialExpression(request, self._adata)

    @responses.activate
    def test_returns_json(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = DifferentialExpression(self.get_request(), self._adata).compute()
        res = res[0].result
        json.loads(res)

    @responses.activate
    def test_returns_a_json_object(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = DifferentialExpression(self.get_request(), self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    @responses.activate
    def test_object_has_all_required_columns(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = DifferentialExpression(self.get_request(), self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for row in res["rows"]:
            keys = sorted(row.keys())
            expected_keys = sorted(
                ["gene_names", "zscore", "abszscore", "qval", "log2fc", "_row"]
            )
            assert keys == expected_keys

    @responses.activate
    def test_appropriate_genes_returned_when_a_limit_is_specified(
        self, mock_dynamo_get
    ):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        request = self.get_request(maxNum=2)

        res = DifferentialExpression(request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)["rows"]

        assert len(res) <= request["body"]["maxNum"]

    @responses.activate
    def test_all_genes_returned_when_no_limit_is_specified(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = DifferentialExpression(self.get_request(), self._adata).compute()
        res = res[0].result
        res = json.loads(res)["rows"]

        assert len(res) <= len(self._adata.raw.var.index)

    @responses.activate
    def test_cells_in_sets_intersection_are_filtered_out(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        MockDynamoClass.setResponse("two_sets_intersected")

        DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="cluster2"), self._adata
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        # Check 1 cell of each of the cell sets is left out
        assert len(baseCells) == len(backgroundCells) == 2

        # Check the cells that haven't been left out are
        # those that are not in the intersection of both sets
        assert len(set(baseCells).intersection(set(backgroundCells))) == 0

    @responses.activate
    def test_cells_not_in_basis_sample_are_filtered_out(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        MockDynamoClass.setResponse("three_sets")

        DifferentialExpression(
            self.get_request(
                cellSet="cluster1", compareWith="cluster2", basis="basisCluster"
            ),
            self._adata,
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        # Check cells not in basis are taken out
        assert len(baseCells) == 1
        assert len(backgroundCells) == 2

    @responses.activate
    def test_rest_keyword_only_adds_cells_in_the_same_hierarchy(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        MockDynamoClass.setResponse("hierarchichal_sets")

        DifferentialExpression(
            self.get_request(cellSet="cluster1", compareWith="rest"),
            self._adata,
        ).compute()

        request_to_r_worker = json.loads(responses.calls[0].request.body)

        baseCells = request_to_r_worker["baseCells"]
        backgroundCells = request_to_r_worker["backgroundCells"]

        # Check there is only one cell in each set
        assert len(baseCells) == 1
        assert len(backgroundCells) == 1
