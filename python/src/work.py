import datetime
from tasks.tasks import TaskFactory
from consume_message import consume
from response import Response
from config import get_config


def main():
    config = get_config()
    last_activity = datetime.datetime.utcnow()
    task_factory = TaskFactory()
    print(datetime.datetime.utcnow(), "Now listening, waiting for work to do...")

    while (
        datetime.datetime.utcnow() - last_activity
    ).total_seconds() <= config.TIMEOUT or config.IGNORE_TIMEOUT:
        msg = consume()
        if msg:
            results = task_factory.submit(msg)
            response = Response(request=msg, results=results)
            response.publish()

            last_activity = datetime.datetime.utcnow()

    print(datetime.datetime.utcnow(), "Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
