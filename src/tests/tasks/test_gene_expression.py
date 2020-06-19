import pytest
import anndata
import os
from tasks.gene_expression import GeneExpression
import json
from botocore.stub import Stubber
import mock
import boto3
from config import get_config
from boto3.dynamodb.types import TypeSerializer

config = get_config()


class TestGeneExpression:
    @pytest.fixture(autouse=True)
    def open_test_adata(self):
        self._adata = anndata.read_h5ad(os.path.join("tests", "test.h5ad"))

    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.correct_request = {
            "experimentId": "5e959f9c9f4b120771249001",
            "timeout": "2099-12-31 00:00:00",
            "body": {
                "name": "GeneExpression",
                "cellSets": ["cluster1"],
                "genes": ["TGFB1", "CST3"],
            },
        }

    """
    Mocks the DynamoDB query for fetching cell sets. Returns a
    correct cell set and yields the patched up object.
    """

    @pytest.fixture
    def mock_dynamo_get(self):
        ser = TypeSerializer()

        response = [
            {
                "name": "my amazing cluster",
                "key": "cluster1",
                "cellIds": ["AAACCGTGCTTCCG-1", "AAAGAGACGCGAGA-1"],
            },
            {
                "name": "my other amazing cluster",
                "key": "cluster2",
                "cellIds": [
                    "TTGGAGACCAATCG-1",
                    "TTGGGAACTGAACC-1",
                    "TTGGTACTCTTAGG-1",
                    "AAAGCAGATATCGG-1",
                ],
            },
        ]

        response = ser.serialize(response)

        test_experiment_id = self.correct_request["experimentId"]
        dynamodb = boto3.resource("dynamodb")
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"cellSets": response}},
            {
                "TableName": config.get_dynamo_table(),
                "Key": {"experimentId": test_experiment_id},
                "ProjectionExpression": "cellSets",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            yield (m, dynamodb)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            GeneExpression()

    def test_throws_on_missing_adata(self):
        with pytest.raises(TypeError):
            GeneExpression(self.correct_request)

    def test_dynamodb_call_is_made_once_when_cell_set_specified(self):
        with mock.patch("boto3.resource") as m:
            global no_called
            no_called = 0

            class MockTable:
                def get_item(*args, **kwargs):
                    global no_called

                    no_called += 1

                    return {"Item": {"cellSets": []}}

            class MockDynamoClass:
                def Table(*args, **kwargs):
                    return MockTable()

            m.return_value = MockDynamoClass()
            GeneExpression(self.correct_request, self._adata).compute()

            assert no_called == 1

    def test_dynamodb_call_is_not_made_when_cell_set_is_all(self):
        self.correct_request["body"]["cellSets"] = "all"

        with mock.patch("boto3.resource") as m:
            global no_called
            no_called = 0

            class MockTable:
                def get_item(*args, **kwargs):
                    global no_called

                    no_called += 1

                    return {"Item": {"cellSets": []}}

            class MockDynamoClass:
                def Table(*args, **kwargs):
                    return MockTable()

            m.return_value = MockDynamoClass()

            print(self.correct_request)
            GeneExpression(self.correct_request, self._adata).compute()

            assert no_called == 0

    def test_works_with_request_and_adata(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb
        GeneExpression(self.correct_request, self._adata)

    def test_returns_json(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        json.loads(res)

    def test_returns_a_json_object(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)
        assert isinstance(res, dict)

    def test_object_returns_all_cells_in_clusters(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res["cells"]) == 2

    def test_object_returns_correct_number_of_cells_with_two_cellsets(
        self, mock_dynamo_get
    ):
        self.correct_request["body"]["cellSets"] = ["cluster1", "cluster2"]
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res["cells"]) == 6

    def test_object_returns_appropriate_number_of_genes(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        assert len(res["data"]) == len(self.correct_request["body"]["genes"])

    def test_object_returns_correct_number_of_expression_data(self, mock_dynamo_get):
        m, dynamodb = mock_dynamo_get
        m.return_value = dynamodb

        res = GeneExpression(self.correct_request, self._adata).compute()
        res = res[0].result
        res = json.loads(res)

        for data in res["data"]:
            assert len(data["expression"]) == 2
