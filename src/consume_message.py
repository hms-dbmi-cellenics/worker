import boto3
from botocore.exceptions import ClientError
import json
from config import get_config
import datetime
import dateutil
import pytz

config = get_config()


def _read_sqs_message():
    sqs = boto3.resource("sqs", **config.BOTO_RESOURCE_KWARGS)

    """
    It is possible that the queue was not created by the time
    the worker launches, because the work queue creation (if needed)
    and the Job spawn are on separate promises and work asyncrhonously.
    This is a performance improvement but it causes the race condition above.

    If this is the case, we just return an empty response
    as if we didn't receive a message in this time frame.
    """
    try:
        queue = sqs.get_queue_by_name(QueueName=config.QUEUE_NAME)
    except ClientError as e:
        if e.response["Error"]["Code"] == "AWS.SimpleQueueService.NonExistentQueue":
            return None
        else:
            raise e

    message = queue.receive_messages(WaitTimeSeconds=20)

    if not message:
        return None

    # Try to parse it as JSON
    try:
        message = message[0]
        print(datetime.datetime.utcnow(), message.body)
        body = json.loads(message.body)
        print(datetime.datetime.utcnow(), "Consumed a message from SQS.")
    except Exception as e:
        print(datetime.datetime.utcnow(), "Exception when loading json: ", e)
        return None
    finally:
        message.delete()

    return body


def consume():
    mssg_body = _read_sqs_message()

    if not mssg_body:
        return None

    timeout = mssg_body["timeout"]
    timeout = dateutil.parser.parse(timeout).astimezone(pytz.utc).replace(tzinfo=None)

    if timeout <= datetime.datetime.utcnow():
        print(
            datetime.datetime.utcnow(),
            "Skipping sending task with uuid",
            mssg_body["uuid"],
            "as its timeout of",
            timeout,
            "has expired...",
        )

        return None

    print(datetime.datetime.utcnow(), mssg_body)
    return mssg_body
