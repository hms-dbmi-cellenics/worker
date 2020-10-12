import os
import types


def get_config():
    kube_env = os.getenv("K8S_ENV")
    cluster_env = os.getenv("CLUSTER_ENV")
    queue_name = os.getenv("WORK_QUEUE")
    timeout = int(os.getenv("WORK_TIMEOUT", default="1200"))

    aws_account_id = os.getenv("AWS_ACCOUNT_ID", default="242905224710")
    aws_region = os.getenv("AWS_DEFAULT_REGION", default="eu-west-1")

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
        QUEUE_NAME=queue_name,
        TIMEOUT=timeout,
        AWS_ACCOUNT_ID=aws_account_id,
        AWS_REGION=aws_region,
        BOTO_RESOURCE_KWARGS={"region_name": aws_region},
        DYNAMO_TABLE=f"experiments-{cluster_env}",
        RESULTS_BUCKET=f"worker-results-{cluster_env}",
        SNS_TOPIC=f"work-results-{cluster_env}",
        R_WORKER_URL="http://localhost:4000",
    )

    if cluster_env == "development" or cluster_env == "test":
        config.QUEUE_NAME = "development-queue.fifo"
        config.AWS_ACCOUNT_ID = "000000000000"
        config.BOTO_RESOURCE_KWARGS["aws_access_key_id"] = "my-key"
        config.BOTO_RESOURCE_KWARGS["aws_secret_access_key"] = "my-secret-key"

    if cluster_env == "development":
        config.BOTO_RESOURCE_KWARGS["endpoint_url"] = "http://host.docker.internal:4566"
        config.R_WORKER_URL = "http://r:4000"

    return config
