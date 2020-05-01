import boto3
import datetime
import os
import io
import json
import anndata
import numpy as np

from tasks.task_map import TASK_MAP

DEFAULT_TIMEOUT = 60 * 20
TOPIC_NAME = "work-results-staging"
QUEUE_NAME = os.getenv("WORK_QUEUE", "test-queue")
TIMEOUT = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))
RESULTS_BUCKET = os.getenv("RESULTS_BUCKET", default="worker-results-staging")
ENVIRONMENT = os.getenv("GITLAB_ENVIRONMENT_NAME", default="local")


def load_file(count_matrix_path):
    if ENVIRONMENT != "local":
        bucket, key = count_matrix_path.split("/", 1)
        try:
            client = boto3.client("s3")
            result = io.BytesIO()
            client.download_fileobj(Bucket=bucket, Key=key, Fileobj=result)
            result.seek(0)

            adata = anndata.read_h5ad(result)
        except Exception as e:
            print("Could not get file from S3", e)
            raise e
    else:
        with open("tgfb1-3-control.h5ad", "rb") as f:
            adata = anndata.read_h5ad(f)

    print("File was loaded.")
    return adata


def run_task(body, adata):
    # Get the contents of the schema.
    task = body["task"]
    details = body["details"]
    task_cls = TASK_MAP[task]

    try:
        result = task_cls(adata).consume(details)
        return result
    except Exception as e:
        # do return this though to the api
        raise e


def upload_result(result, bucket, key):
    client = boto3.client("s3")
    r = client.put_object(Key=key, Bucket=bucket, Body=json.dumps(result))
    return r


def publish_message(message):
    sns = boto3.client("sns")
    resp = sns.publish(
        TargetArn="arn:aws:sns:eu-west-2:242905224710:work-results-staging",
        Message=json.dumps({"default": json.dumps(message)}),
        Subject="a short subject for your message",
        MessageStructure="json",
    )
    print("Result was published.")
    print(resp)


def main():
    """
    Consumes and parses instructions from a message queue.
    """
    sqs = boto3.resource("sqs")
    queue = sqs.get_queue_by_name(QueueName=QUEUE_NAME)
    last_activity = datetime.datetime.now()
    adata = None

    print("Now listening, waiting for work to do...")

    while (datetime.datetime.now() - last_activity).total_seconds() <= TIMEOUT:
        message = queue.receive_messages(WaitTimeSeconds=20)

        if not message:
            continue

        # Try to parse it as JSON
        try:
            message = message[0]
            body = json.loads(message.body)
            print("Consumed a message from S3.")
        except Exception as e:
            print("Exception when loading json: ", e)
            continue
        finally:
            message.delete()

        # Load file from S3 if not already present.
        if not adata:
            count_matrix_path = body["count_matrix"]
            adata = load_file(count_matrix_path)

        # Run task.
        result = run_task(body, adata)
        result_json = {"result": result.tolist()}

        uuid = body.get("uuid", "1234")
        s3_path = "/".join([RESULTS_BUCKET, uuid])

        resp_body = {
            "uuid": uuid,
            "socketId": body.get("socketId", "567"),
            "results": [
                {
                    "content-type": "application/json",
                    "type": "s3-path",
                    "body": s3_path,
                }
            ],
        }

        print(resp_body)

        # Publish task results
        upload_result(result_json, RESULTS_BUCKET, uuid)
        publish_message(resp_body)

        last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
