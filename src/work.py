import boto3
import datetime
import os
import io
import json
import anndata

from tasks.task_map import TASK_MAP

DEFAULT_TIMEOUT = 60 * 20
TOPIC_NAME = "analysis-results"
QUEUE_NAME = os.getenv("WORK_QUEUE", "test-queue")
TIMEOUT = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))


def load_file(count_matrix_path):
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


def main():
    """
    Consumes and parses instructions from a message queue.

    It takes the following environment variables: WORK_QUEUE (the SQS queue
    to join), WORK_TIMEOUT (the timeout for this session in
    seconds, default: 60*20).
    """
    queue = os.getenv("WORK_QUEUE", "test-queue")

    if not queue:
        raise ValueError("No work queue specified.")

    sqs = boto3.resource("sqs")
    queue = sqs.get_queue_by_name(QueueName=QUEUE_NAME)
    print(f"Now listening, waiting for work to do...")

    last_activity = datetime.datetime.now()
    adata = None

    while (datetime.datetime.now() - last_activity).total_seconds() <= TIMEOUT:
        message = queue.receive_messages(WaitTimeSeconds=20)

        if not message:
            continue

        message = message[0]

        # Try to parse it as JSON
        try:
            body = json.loads(message.body)
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

        # Publish task results to SNS topic
        print(result)

        last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
