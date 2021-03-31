import datetime
from tasks.tasks import TaskFactory
from consume_message import consume
from response import Response
from config import get_config
from aws_xray_sdk.core import xray_recorder
from aws_xray_sdk import global_sdk_config


def main():
    global_sdk_config.set_sdk_enabled(False)

    config = get_config()
    last_activity = datetime.datetime.utcnow()
    task_factory = TaskFactory()
    print(datetime.datetime.utcnow(), "Now listening, waiting for work to do...")

    if config.IGNORE_TIMEOUT:
        print(
            datetime.datetime.utcnow(),
            "Worker configured to ignore timeout, will run forever...",
        )

    while (
        datetime.datetime.utcnow() - last_activity
    ).total_seconds() <= config.TIMEOUT or config.IGNORE_TIMEOUT:
        global_sdk_config.set_sdk_enabled(False)

        msg = consume()
        if msg:
            results = task_factory.submit(msg)
            response = Response(request=msg, results=results)
            response.publish()

            last_activity = datetime.datetime.utcnow()

        xray_recorder.end_segment()

    print(datetime.datetime.utcnow(), "Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
