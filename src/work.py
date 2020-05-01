import boto3
import datetime
import os
import io
import json
import anndata

from tasks.task_map import TASK_MAP
from consume_message import consume
from result import Result

DEFAULT_TIMEOUT = 60 * 20
TOPIC_NAME = "work-results-staging"
TIMEOUT = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))


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
    last_activity = datetime.datetime.now()
    adata = None
    print("Now listening, waiting for work to do...")

    while (datetime.datetime.now() - last_activity).total_seconds() <= TIMEOUT:
        adata, body = consume(adata)
        r = run_task(body, adata)
        result = Result(work_def=body, result=r)
        result.publish()

        last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
