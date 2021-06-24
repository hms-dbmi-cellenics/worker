import boto3
from botocore.stub import Stubber, ANY
import mock
from worker.config import config
from worker.consume_message import _read_sqs_message, consume


class TestConsumeMessage:
    def test_read_sqs_message_fetches_messages_from_the_correct_queue(self):
        sqs = boto3.resource("sqs", **config.BOTO_RESOURCE_KWARGS)
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
                "AttributeNames": ["AWSTraceHeader"],
            },
        )

        with mock.patch("boto3.resource") as m, stubber:
            m.return_value = sqs
            _read_sqs_message()

    def test_read_sqs_message_returns_falsy_on_non_existent_queue(self):
        sqs = boto3.resource("sqs", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(sqs.meta.client)
        stubber.add_client_error(
            "get_queue_url",
            service_error_code="AWS.SimpleQueueService.NonExistentQueue",
            http_status_code=400,
            expected_params={"QueueName": config.QUEUE_NAME},
        )

        with mock.patch("boto3.resource") as m, stubber:
            m.return_value = sqs
            r = consume()
            assert not r

    def test_read_sqs_message_returns_falsy_on_no_incoming_message(self):
        sqs = boto3.resource("sqs", **config.BOTO_RESOURCE_KWARGS)
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
                "AttributeNames": ["AWSTraceHeader"],
            },
        )

        with mock.patch("boto3.resource") as m, stubber:
            m.return_value = sqs
            r = _read_sqs_message()

            assert not r

    def test_read_sqs_message_returns_falsy_on_badly_formatted_message(self):
        sqs = boto3.resource("sqs", **config.BOTO_RESOURCE_KWARGS)
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
                    {"MessageId": "asd", "ReceiptHandle": "ewrwe", "Body": '{"not_json_asd: "b"}'}
                ]
            },
            {
                "QueueUrl": "my_very_valid_and_existing_queue_url",
                "WaitTimeSeconds": ANY,
                "AttributeNames": ["AWSTraceHeader"],
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

        with mock.patch("boto3.resource") as m, stubber:
            m.return_value = sqs
            r = _read_sqs_message()

            assert not r

    def test_request_with_expired_timeout_is_discarded(self):
        request = {
            "experimentId": "random-experiment-id",
            "timeout": "2000-01-01 00:00:00",
            "uuid": "random-uuid",
        }

        with mock.patch("worker.consume_message._read_sqs_message") as m:
            m.return_value = request
            result = consume()

            assert result is None

    def test_consume_request_with_non_expired_timeout_successfully(self):
        request = {
            "experimentId": "random-experiment-id",
            "timeout": "2900-01-01 00:00:00",
            "uuid": "random-uuid",
        }

        with mock.patch("worker.consume_message._read_sqs_message") as m:
            m.return_value = request
            result = consume()

            assert result == {
                "experimentId": "random-experiment-id",
                "timeout": "2900-01-01 00:00:00",
                "uuid": "random-uuid",
            }
