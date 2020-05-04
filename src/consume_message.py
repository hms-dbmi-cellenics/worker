import boto3
import os
import io
import json
import anndata

QUEUE_NAME = os.getenv("WORK_QUEUE", "test-queue")
ENVIRONMENT = os.getenv("GITLAB_ENVIRONMENT_NAME", default="local")
DYNAMO_TABLE = os.getenv("DYNAMODB_TABLE", default="experiments-staging")


def _read_sqs_message():
    sqs = boto3.resource("sqs")
    queue = sqs.get_queue_by_name(QueueName=QUEUE_NAME)
    message = queue.receive_messages(WaitTimeSeconds=20)

    if not message:
        return None

    # Try to parse it as JSON
    try:
        message = message[0]
        print(message.body)
        body = json.loads(message.body)
        print("Consumed a message from S3.")
    except Exception as e:
        print("Exception when loading json: ", e)
        return None
    finally:
        message.delete()

    return body


def _get_matrix_path(experiment_id):
    dynamo = boto3.resource("dynamodb").Table(DYNAMO_TABLE)
    matrix_path = dynamo.get_item(
        Key={"experimentId": experiment_id}, ProjectionExpression="matrixPath",
    )
    print("successfully got the matrix path from database.")
    return matrix_path


def _load_file(matrix_path):
    if ENVIRONMENT != "local":
        bucket, key = matrix_path.split("/", 1)
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


def consume(adata):
    mssg_body = _read_sqs_message()

    if not mssg_body:
        return adata, None

    if not adata:
        print("adata does not exist, I need to download it ...")
        matrix_path = _get_matrix_path(mssg_body["experimentId"])
        adata = _load_file(matrix_path)

    print(mssg_body)
    return adata, mssg_body
