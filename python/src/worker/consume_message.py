import datetime
import json
import traceback
from logging import info

import aws_xray_sdk as xray
import boto3
import dateutil
import pytz
from aws_xray_sdk.core import xray_recorder
from aws_xray_sdk.core.models.trace_header import TraceHeader
from botocore.exceptions import ClientError

from .config import config


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

    message = queue.receive_messages(
        WaitTimeSeconds=20, AttributeNames=["AWSTraceHeader"]
    )

    if not message:
        return None

    # Try to parse it as JSON
    try:
        message = message[0]

        trace_header = message.attributes and message.attributes.get(
            "AWSTraceHeader", None
        )

        if trace_header:
            xray.global_sdk_config.set_sdk_enabled(True)

            header = TraceHeader.from_header_str(trace_header)
            trace_id = header.root
            sampled = header.sampled

            xray_recorder.begin_segment(
                f"worker-{config.CLUSTER_ENV}-{config.SANDBOX_ID}",
                traceid=trace_id,
                sampling=sampled,
                parent_id=header.parent,
            )

        body = json.loads(message.body)
        info("Consumed a message from SQS.")
    except Exception as e:
        xray_recorder.current_segment().add_exception(e, traceback.format_exc())

        info("Exception when loading message", message.body)
        info("Exception:", e)
        return None
    finally:
        message.delete()

    return body


@xray_recorder.capture("consume_message._response_exists")
def _response_exists(mssg_body):
    client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
    ETag = mssg_body["ETag"]

    print("&&&&&&&&&&& ", client)
    response = client.list_objects_v2(
        Bucket=config.RESULTS_BUCKET,
        Prefix=ETag,
    )

    print("^^^^^^^^^^^^ ", response)
    for obj in response.get("Contents", []):
        print("++++++ ", obj)
        if obj["Key"] == ETag:
            return obj["Size"]


def consume():
    mssg_body = _read_sqs_message()

    if not mssg_body:
        return None

    timeout = mssg_body["timeout"]
    timeout = dateutil.parser.parse(timeout).astimezone(pytz.utc).replace(tzinfo=None)

    if timeout <= datetime.datetime.utcnow():
        info(
            f"Skipping processing task with ETag {mssg_body['ETag']} "
            f"as its timeout of {timeout} has expired..."
        )

        return None

    if _response_exists(mssg_body):
        info(
            f"Skipping processing task with ETag {mssg_body['ETag']} "
            f"as a response with this hash is already in S3."
        )
        print("IT exists")
        return None

    print("IT doesn't exist")
    info(json.dumps(mssg_body, indent=2, sort_keys=True))
    return mssg_body
