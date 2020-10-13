import datetime
from tasks.tasks import TaskFactory
from consume_message import consume
from response import Response
from config import get_config
from helpers.r_readiness import check_r_readiness


def main():
    config = get_config()

    check_r_readiness()

    tf = TaskFactory()

    last_activity = datetime.datetime.utcnow()
    adata = None
    print(datetime.datetime.utcnow(), "Now listening, waiting for work to do...")

    while (
        datetime.datetime.utcnow() - last_activity
    ).total_seconds() <= config.TIMEOUT:
        msg = consume()
        if msg:
            results, adata = tf.submit(msg, adata)
            if not adata:
                raise Exception("Adata file did not get loaded properly")
            response = Response(request=msg, results=results)
            response.publish()

            last_activity = datetime.datetime.utcnow()

    print(datetime.datetime.utcnow(), "Timeout exceeded, shutting down...")


if __name__ == "__main__":
    main()
