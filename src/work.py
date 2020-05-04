import datetime
from tasks.tasks import TaskFactory
from consume_message import consume
from result import Result
from config import get_config


def main():
    config = get_config()
    last_activity = datetime.datetime.now()
    adata = None
    print("Now listening, waiting for work to do...")

    while (datetime.datetime.now() - last_activity).total_seconds() <= config.TIMEOUT:
        adata, mssg = consume(adata)
        if mssg:
            r = TaskFactory().submit(mssg["body"], adata)
            result = Result(work_def=mssg, result=r)
            result.publish()
            last_activity = datetime.datetime.now()

    print("Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
