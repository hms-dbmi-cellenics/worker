import os
import types


def get_config():
    gitlab_env = os.getenv("GITLAB_ENVIRONMENT_NAME")
    cluster_env = os.getenv("CLUSTER_ENV")
    queue_name = os.getenv("WORK_QUEUE")
    timeout = int(os.getenv("WORK_TIMEOUT", default="1200"))

    aws_account_id = os.getenv("AWS_ACCOUNT_ID", default="242905224710")
    aws_region = os.getenv("AWS_DEFAULT_REGION", default="eu-west-2")

    # set up cluster env based on gitlab env if one was not specified
    # this is only run if `gitlab_env` is specified, i.e. when the system
    # is run in staging/production or in testing
    if gitlab_env and not cluster_env:
        if gitlab_env == "production":
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
    )

    if cluster_env == "development":
        config.QUEUE_NAME = "development-queue.fifo"
        config.AWS_ACCOUNT_ID = "000000000000"
        config.BOTO_RESOURCE_KWARGS["endpoint_url"] = "http://localhost:4566"

    return config
