import boto3
import datetime
import os
import io
import json
import anndata

from tasks.task_map import TASK_MAP

DEFAULT_TIMEOUT = 60 * 20


def load_file(bucket, key):
    client = boto3.client("s3")

    result = io.BytesIO()
    client.download_fileobj(Bucket=bucket, Key=key, Fileobj=result)
    result.seek(0)

    adata = anndata.read_h5ad(result)

    return adata


def main():
    """
    Consumes and parses instructions from a message queue.

    It takes the following environment variables: WORK_QUEUE (the SQS queue
    to join), WORK_TIMEOUT (the timeout for this session in
    seconds, default: 60*20).
    """

    queue = os.getenv("WORK_QUEUE", None)

    if not queue:
        raise ValueError("No work queue specified.")

    try:
        timeout = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))
    except Exception:
        raise ValueError(
            "WORK_TIMEOUT is set to an invalid value (must be an integer)."
        )

    sqs = boto3.resource("sqs")
    queue = sqs.get_queue_by_name(QueueName=queue)
    print(f"Now listening, waiting for work to do...")

    last_activity = datetime.datetime.now()
    adata = None

    while (datetime.datetime.now() - last_activity).total_seconds() <= timeout:
        message = queue.receive_messages(WaitTimeSeconds=20)

        if not message:
            continue

        message = message[0]

        # Try to parse it as JSON
        try:
            body = json.loads(message.body)
        except Exception:
            print("Malformed JSON:", message.body)
            continue
        finally:
            message.delete()

        # Get the contents of the schema.
        try:
            count_matrix_path = body["count_matrix"]
            task = body["task"]
            details = body["details"]
        except Exception as e:
            print("Illegal request", e)

        # Load file from S3 if not already present.
        if not adata:
            try:
                s3_path = count_matrix_path.split("/", 1)
                adata = load_file(*s3_path)
            except Exception as e:
                print("Could not get file from S3", e)
                raise e

        # Find the task function from the map.
        try:
            task_cls = TASK_MAP[task]
        except Exception as e:
            print("Task {} not recognized:".format(task), e)

        # Run task.
        try:
            result = task_cls(adata).consume(details)
        except Exception as e:
            # do return this though to the api
            raise e

        print(result)

        last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
