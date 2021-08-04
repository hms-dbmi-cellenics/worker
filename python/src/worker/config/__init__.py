import os
import types
import re

from aws_xray_sdk import core, global_sdk_config


def get_label(label_key):
    labels = {}

    try:
        with open("labels") as f:
            for line in f.readlines():
                key, value = line.rstrip("\n").replace('"', "").split("=")
                labels[key] = value
    except FileNotFoundError:
        pass
    
    # Attempt to get the data directly from the label. If the label
    # does not exist (because e.g. it is in development or because
    # the worker is unassigned to an experiment) we try to get the
    # info from an env variable (experimentId -> EXPERIMENT_ID).
    # If unsuccessful, we return None.
    return labels.get(
        label_key,
        os.getenv(
            re.sub(r'(?<!^)(?=[A-Z])', '_', label_key).upper(),
            None
        )
    )

kube_env = os.getenv("K8S_ENV")
cluster_env = os.getenv("CLUSTER_ENV")
queue_name = os.getenv("WORK_QUEUE")

# timeout is in seconds, set to 1 hour
timeout = int(os.getenv("WORK_TIMEOUT", default=str(60 * 60 * 9)))
ignore_timeout = os.getenv("IGNORE_TIMEOUT") == "true"

aws_account_id = os.getenv("AWS_ACCOUNT_ID", default="242905224710")
aws_region = os.getenv("AWS_DEFAULT_REGION", default="eu-west-1")

# set up cluster env based on github env if one was not specified
# this is only run if `kube_env` is specified, i.e. when the system
# is run in staging/production or in testing
if kube_env and not cluster_env:
    if kube_env == "production":
        cluster_env = "production"
    else:
        cluster_env = "staging"

if not cluster_env:
    cluster_env = "development"

config = types.SimpleNamespace(
    CLUSTER_ENV=cluster_env,
    SANDBOX_ID=property(lambda: get_label('sandboxId'), lambda: None),
    QUEUE_NAME=queue_name,
    TIMEOUT=timeout,
    IGNORE_TIMEOUT=ignore_timeout,
    AWS_ACCOUNT_ID=aws_account_id,
    AWS_REGION=aws_region,
    BOTO_RESOURCE_KWARGS={"region_name": aws_region},
    DYNAMO_TABLE=f"experiments-{cluster_env}",
    CELL_SETS_BUCKET=f"cell-sets-{cluster_env}",
    SOURCE_BUCKET=f"processed-matrix-{cluster_env}",
    RESULTS_BUCKET=f"worker-results-{cluster_env}",
    SNS_TOPIC=property(lambda: f"work-results-{cluster_env}-{get_label('experimentId')}", lambda: None),
    R_WORKER_URL="http://localhost:4000",
    EXPERIMENT_ID=property(lambda: get_label('experimentId'), lambda: None),
    # this works because in CI, `data/` is deployed under `worker/`
    # whereas in a container, it is mounted to `/data`. Either way, this ensures
    # that the appropriate path is selected, as both are two directories up
    LOCAL_DIR=os.path.join(os.pardir, os.pardir, "data"),
)

if cluster_env == "development" or cluster_env == "test":
    config.QUEUE_NAME = "development-queue.fifo"
    config.AWS_ACCOUNT_ID = "000000000000"
    config.BOTO_RESOURCE_KWARGS["aws_access_key_id"] = "my-key"
    config.BOTO_RESOURCE_KWARGS["aws_secret_access_key"] = "my-secret-key"
    global_sdk_config.set_sdk_enabled(False)

if cluster_env == "development":
    config.BOTO_RESOURCE_KWARGS[
        "endpoint_url"
    ] = "http://host.docker.internal:4566"
    config.R_WORKER_URL = "http://r:4000"

if cluster_env != "test":
    core.patch_all()
