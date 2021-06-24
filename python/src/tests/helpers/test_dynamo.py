import boto3
import mock
from botocore.stub import ANY, Stubber

from worker.config import config
from worker.helpers.dynamo import get_item_from_dynamo


class TestDynamo:
    def test_get_matrix_path_attempts_connection_to_appropriate_table(self):
        dynamodb = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS)

        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.DYNAMO_TABLE,
                "Key": ANY,
                "ProjectionExpression": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            resp = get_item_from_dynamo("my-very-serious-experiment", "matrixPath")
            assert resp == "very/genuine/path"

    def test_get_matrix_path_gets_correct_experiment_id_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.DYNAMO_TABLE,
                "Key": {"experimentId": test_experiment_id},
                "ProjectionExpression": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            resp = get_item_from_dynamo(test_experiment_id, "matrixPath")
            assert resp == "very/genuine/path"

    def test_get_matrix_path_gets_correct_field_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.DYNAMO_TABLE,
                "Key": ANY,
                "ProjectionExpression": "matrixPath",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            resp = get_item_from_dynamo(test_experiment_id, "matrixPath")
            assert resp == "very/genuine/path"

    def test_get_non_existing_field_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {}},
            {
                "TableName": config.DYNAMO_TABLE,
                "Key": ANY,
                "ProjectionExpression": "matrixPath",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            resp = get_item_from_dynamo(test_experiment_id, "matrixPath")
            assert resp == {}

    def test_get_non_existing_experiment_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {},
            {
                "TableName": config.DYNAMO_TABLE,
                "Key": ANY,
                "ProjectionExpression": "matrixPath",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            resp = get_item_from_dynamo(test_experiment_id, "matrixPath")
            assert resp == {}
