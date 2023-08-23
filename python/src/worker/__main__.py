import datetime
import time
from logging import INFO, basicConfig, info

import aws_xray_sdk as xray
from aws_xray_sdk.core import xray_recorder

from .config import config
from .consume_message import consume
from .response import Response
from .tasks.factory import TaskFactory
from socket_io_emitter import Emitter
from .constants import STARTED_TASK

# configure logging
basicConfig(format="%(asctime)s %(message)s", level=INFO)


def main():
    if config.IGNORE_TIMEOUT:
        info("Worker configured to ignore timeout, will run forever...")

    while not config.EXPERIMENT_ID:
        info("Experiment not yet assigned, waiting...")
        time.sleep(5)

    # Disable X-Ray for initial setup so we don't end up
    # with segment warnings before any message is sent
    xray.global_sdk_config.set_sdk_enabled(False)

    last_activity = datetime.datetime.utcnow()
    task_factory = TaskFactory()
    info(
        f"Now listening for experiment {config.EXPERIMENT_ID}, waiting for work to do..."
    )

    while (
        datetime.datetime.utcnow() - last_activity
    ).total_seconds() <= config.TIMEOUT or config.IGNORE_TIMEOUT:

        # Disable X-Ray before message is identified and processed
        xray.global_sdk_config.set_sdk_enabled(False)

        request = consume()
        if request:
            io = Emitter({"client": config.REDIS_CLIENT})
            io.Emit(f'Heartbeat-{request["experimentId"]}', {"type": "WorkResponse", "workingOn": STARTED_TASK, "request": request})

            result = task_factory.submit(request)

            response = Response(request, result)
            response.publish()

            last_activity = datetime.datetime.utcnow()

        xray_recorder.end_segment()

    info("Timeout exceeded, shutting down...")


main()
