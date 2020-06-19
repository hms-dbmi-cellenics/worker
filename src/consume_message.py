import boto3
from botocore.exceptions import ClientError
import io
import json
import anndata
from config import get_config
import datetime
import dateutil
import pytz

config = get_config()


def _read_sqs_message():
    sqs = boto3.resource("sqs", region_name=config.AWS_REGION)

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
        print(datetime.datetime.now(), message.body)
        body = json.loads(message.body)
        print(datetime.datetime.now(), "Consumed a message from S3.")
    except Exception as e:
        print(datetime.datetime.now(), "Exception when loading json: ", e)
        return None
    finally:
        message.delete()

    return body


def _get_matrix_path(experiment_id):
    dynamo = boto3.resource("dynamodb", region_name=config.AWS_REGION).Table(
        config.get_dynamo_table()
    )

    # todo: the projectionexpression stopped working for some reason, fix it!
    resp = dynamo.get_item(
        Key={"experimentId": experiment_id}, ProjectionExpression="matrixPath",
    )
    matrix_path = resp["Item"]["matrixPath"]
    print("successfully got the matrix path from database.")
    return matrix_path


def _load_file(matrix_path):
    print(config.ENVIRONMENT)
    if config.ENVIRONMENT != "development":
        print(datetime.datetime.now(), "Have to download anndata file from s3")
        bucket, key = matrix_path.split("/", 1)
        try:
            client = boto3.client("s3")
            result = io.BytesIO()
            client.download_fileobj(Bucket=bucket, Key=key, Fileobj=result)
            result.seek(0)

            adata = anndata.read_h5ad(result)
        except Exception as e:
            print(datetime.datetime.now(), "Could not get file from S3", e)
            raise e
    else:
        with open("./tests/test.h5ad", "rb") as f:
            adata = anndata.read_h5ad(f)

    print(datetime.datetime.now(), "File was loaded.")
    return adata


def consume(adata):
    mssg_body = _read_sqs_message()

    if not mssg_body:
        return adata, None

    timeout = mssg_body["timeout"]
    timeout = dateutil.parser.parse(timeout).astimezone(pytz.utc).replace(tzinfo=None)

    if timeout <= datetime.datetime.now():
        print(
            datetime.datetime.now(),
            "Skipping sending task with uuid",
            mssg_body["uuid"],
            "as its timeout of",
            timeout,
            "has expired...",
        )

        return adata, None

    if not adata:
        print(
            datetime.datetime.now(), "adata does not exist, I need to download it ..."
        )
        matrix_path = _get_matrix_path(mssg_body["experimentId"])
        adata = _load_file(matrix_path)

    print(datetime.datetime.now(), mssg_body)
    return adata, mssg_body
