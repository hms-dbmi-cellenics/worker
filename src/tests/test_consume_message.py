import boto3
from botocore.stub import Stubber, ANY
import mock
from config import get_config
from consume_message import _read_sqs_message, _get_matrix_path, _load_file, consume
from moto import mock_s3

config = get_config()


class TestConsumeMessage:
    def test_read_sqs_message_fetches_messages_from_the_correct_queue(self):
        sqs = boto3.resource("sqs")
        stubber = Stubber(sqs.meta.client)
        stubber.add_response(
            "get_queue_url",
            {"QueueUrl": "my_very_valid_and_existing_queue_url"},
            {"QueueName": config.QUEUE_NAME},
        )
        stubber.add_response(
            "receive_message",
            {},
            {
                "QueueUrl": "my_very_valid_and_existing_queue_url",
                "WaitTimeSeconds": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = sqs
            _read_sqs_message()

    def test_read_sqs_message_returns_falsy_on_non_existent_queue(self):
        sqs = boto3.resource("sqs")
        stubber = Stubber(sqs.meta.client)
        stubber.add_client_error(
            "get_queue_url",
            service_error_code="AWS.SimpleQueueService.NonExistentQueue",
            http_status_code=400,
            expected_params={"QueueName": config.QUEUE_NAME},
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = sqs
            r = _read_sqs_message()

            assert not r

    def test_read_sqs_message_returns_falsy_on_no_incoming_message(self):
        sqs = boto3.resource("sqs")
        stubber = Stubber(sqs.meta.client)
        stubber.add_response(
            "get_queue_url",
            {"QueueUrl": "my_very_valid_and_existing_queue_url"},
            {"QueueName": config.QUEUE_NAME},
        )
        stubber.add_response(
            "receive_message",
            {"Messages": []},
            {
                "QueueUrl": "my_very_valid_and_existing_queue_url",
                "WaitTimeSeconds": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = sqs
            r = _read_sqs_message()

            assert not r

    def test_read_sqs_message_returns_falsy_on_badly_formatted_message(self):
        sqs = boto3.resource("sqs")
        stubber = Stubber(sqs.meta.client)
        stubber.add_response(
            "get_queue_url",
            {"QueueUrl": "my_very_valid_and_existing_queue_url"},
            {"QueueName": config.QUEUE_NAME},
        )
        stubber.add_response(
            "receive_message",
            {
                "Messages": [
                    {"MessageId": "asd", "ReceiptHandle": "ewrwe", "Body": '{"a": "b"}'}
                ]
            },
            {
                "QueueUrl": "my_very_valid_and_existing_queue_url",
                "WaitTimeSeconds": ANY,
            },
        )
        stubber.add_response(
            "delete_message",
            {},
            {
                "QueueUrl": "my_very_valid_and_existing_queue_url",
                "ReceiptHandle": "ewrwe",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = sqs
            r = _read_sqs_message()

            assert r == {"a": "b"}

    def test_get_matrix_path_attempts_connection_to_appropriate_table(self):
        dynamodb = boto3.resource("dynamodb")
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.get_dynamo_table(),
                "Key": ANY,
                "ProjectionExpression": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            _get_matrix_path("my-very-serious-experiment")

    def test_get_matrix_path_gets_correct_experiment_id_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb")
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.get_dynamo_table(),
                "Key": {"experimentId": test_experiment_id},
                "ProjectionExpression": ANY,
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            _get_matrix_path(test_experiment_id)

    def test_get_matrix_path_gets_correct_field_from_db(self):
        test_experiment_id = "my-very-serious-experiment"
        dynamodb = boto3.resource("dynamodb")
        stubber = Stubber(dynamodb.meta.client)
        stubber.add_response(
            "get_item",
            {"Item": {"matrixPath": {"S": "very/genuine/path"}}},
            {
                "TableName": config.get_dynamo_table(),
                "Key": ANY,
                "ProjectionExpression": "matrixPath",
            },
        )
        stubber.activate()

        with mock.patch("boto3.resource") as m:
            m.return_value = dynamodb
            _get_matrix_path(test_experiment_id)

    @mock_s3
    def test_load_file_returns_correct_adata_object_when_path_and_key_exist_in_non_development(
        self,
    ):
        s3 = boto3.client("s3")
        bucket = "my_custom_bucket_path"
        key = "very/long/and/convoluted/path"
        s3.create_bucket(Bucket=bucket)

        with open("tests/test.h5ad", "rb") as f, mock.patch("config.get_config") as m:
            mock_config = config
            mock_config.ENVIRONMENT = "staging"
            m.return_value = mock_config

            s3.upload_fileobj(f, bucket, key)
            a = _load_file(f"{bucket}/{key}")

            assert "AnnData" in type(a).__name__

    def test_request_with_expired_timeout_is_discarded(self):
        request = {
            "experimentId": "random-experiment-id",
            "timeout": "2000-01-01 00:00:00",
            "uuid": "random-uuid",
        }

        with open("tests/test.h5ad", "rb") as f, mock.patch(
            "consume_message._read_sqs_message"
        ) as m:
            m.return_value = request
            result = consume(f)

            assert result == (f, None)
