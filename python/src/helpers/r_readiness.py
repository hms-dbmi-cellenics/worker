import datetime
import requests
import backoff
from config import get_config

config = get_config()


@backoff.on_exception(
    backoff.expo,
    Exception,
    max_time=20,
    giveup=lambda e: e.response and e.response.status_code >= 400,
)
def check_r_readiness():
    print(datetime.datetime.utcnow(), "Attempting to connect to R worker...")

    r = requests.head(f"{config.R_WORKER_URL}/health")
    r.raise_for_status()

    print(
        datetime.datetime.utcnow(),
        f"Successful connection, status code {r.status_code}.",
    )
