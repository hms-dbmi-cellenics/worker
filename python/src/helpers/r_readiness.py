import datetime
import requests
import backoff
from config import get_config

config = get_config()


def fatal_error(e):
    if e.response and 400 <= e.response.status_code:
        return True


@backoff.on_exception(backoff.expo, Exception, max_time=20, giveup=fatal_error)
def check_r_readiness():
    print(datetime.datetime.utcnow(), "Attempting to connect to R worker...")

    r = requests.head(f"{config.R_WORKER_URL}/health")
    r.raise_for_status()
