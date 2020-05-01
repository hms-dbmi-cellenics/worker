import datetime
import os
from tasks.tasks import TaskFactory
from consume_message import consume
from result import Result

DEFAULT_TIMEOUT = 60 * 20
TOPIC_NAME = "work-results-staging"
TIMEOUT = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))


def run_task(adata, body):
    # Get the contents of the schema.
    print("++++", body)
    task_type = body["task"]
    details = body["details"]

    try:
        result = TaskFactory.factory(task_type, adata).compute(details)
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
        r = run_task(adata, body)
        result = Result(work_def=body, result=r)
        result.publish()
        last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
