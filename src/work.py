import datetime
import os
from tasks.tasks import TaskFactory
from consume_message import consume
from result import Result

DEFAULT_TIMEOUT = 60 * 20
TIMEOUT = int(os.getenv("WORK_TIMEOUT", default=DEFAULT_TIMEOUT))


def main():
    last_activity = datetime.datetime.now()
    adata = None
    print("Now listening, waiting for work to do...")

    while (datetime.datetime.now() - last_activity).total_seconds() <= TIMEOUT:
        adata, mssg = consume(adata)
        if mssg:
            body = mssg["body"]
            r = TaskFactory().submit(body, adata)
            result = Result(work_def=body, result=r)
            result.publish()
            last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
