import datetime
from tasks.tasks import TaskFactory
from consume_message import consume
from response import Response
from config import get_config


def main():
    config = get_config()
    last_activity = datetime.datetime.now()
    adata = None
    print(datetime.datetime.now(), "Now listening, waiting for work to do...")

    while (datetime.datetime.now() - last_activity).total_seconds() <= config.TIMEOUT:
        msg = consume()
        if msg:
            results = TaskFactory().submit(msg, adata)

            response = Response(request=msg, results=results)
            response.publish()

            last_activity = datetime.datetime.now()

    print(datetime.datetime.now(), "Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
