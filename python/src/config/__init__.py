import os
import types
import aws_xray_sdk as xray

def get_config():
    kube_env = os.getenv("K8S_ENV")
    cluster_env = os.getenv("CLUSTER_ENV")
    queue_name = os.getenv("WORK_QUEUE")
    sandbox_id = os.getenv("SANDBOX_ID", default="default")

    # timeout is in seconds, set to 1 hour
    timeout = int(os.getenv("WORK_TIMEOUT", default="3600"))
    ignore_timeout = os.getenv("IGNORE_TIMEOUT") == "true"

    aws_account_id = os.getenv("AWS_ACCOUNT_ID", default="242905224710")
    aws_region = os.getenv("AWS_DEFAULT_REGION", default="eu-west-1")
    experiment_id = os.getenv(
        "EXPERIMENT_ID", default="e52b39624588791a7889e39c617f669e"
    )

    # set up cluster env based on gitlab env if one was not specified
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
        SANDBOX_ID=sandbox_id,
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
        SNS_TOPIC=f"work-results-{cluster_env}-{sandbox_id}",
        R_WORKER_URL="http://localhost:4000",
        EXPERIMENT_ID=experiment_id,
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
        xray.global_sdk_config.set_sdk_enabled(False)

    if cluster_env == "development":
        config.BOTO_RESOURCE_KWARGS["endpoint_url"] = "http://host.docker.internal:4566"
        config.R_WORKER_URL = "http://r:4000"

    if cluster_env != "test":
        xray.core.patch_all()

    return config
