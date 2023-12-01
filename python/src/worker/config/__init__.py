import os
import re
import types
from functools import cached_property
from .domain_specific import get_domain_specific

import boto3
import redis
from aws_xray_sdk import core, global_sdk_config


kube_env = os.getenv("K8S_ENV")
cluster_env = os.getenv("CLUSTER_ENV")
timeout = get_domain_specific().get(kube_env, {}).get('timeout', 10 * 60)

ignore_timeout = os.getenv("IGNORE_TIMEOUT") == "true"
aws_account_id = os.getenv("AWS_ACCOUNT_ID")
aws_region = os.getenv("AWS_DEFAULT_REGION")

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


class Config(types.SimpleNamespace):
    @staticmethod
    def get_label(label_key, default=None):
        labels = {}

        try:
            with open("/etc/podinfo/labels") as f:
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
            os.getenv(re.sub(r"(?<!^)(?=[A-Z])", "_", label_key).upper(), default),
        )

    @property
    def EXPERIMENT_ID(self):
        label = self.get_label("experimentId")
        if label is None and cluster_env == "development":
            return "test-experiment-id"
        return label

    @property
    def SANDBOX_ID(self):
        return self.get_label("sandboxId")

    @property
    def WORK_QUEUE_HASH(self):
        return self.get_label("workQueueHash")

    @property
    def QUEUE_NAME(self):
        if cluster_env == "development" or cluster_env == "test":
            return "development-queue.fifo"

        return f"queue-job-{self.WORK_QUEUE_HASH}-{cluster_env}.fifo"

    @cached_property
    def REDIS_ENDPOINT(self):
        client = boto3.client("elasticache", **self.BOTO_RESOURCE_KWARGS)
        response = client.describe_replication_groups(
            ReplicationGroupId=f"biomage-redis-{cluster_env}"
        )

        # only one group matches the ID
        replication_group = response["ReplicationGroups"][0]

        # only one node group for redis
        node_group = replication_group["NodeGroups"][0]

        return node_group["PrimaryEndpoint"]

    @cached_property
    def REDIS_CLIENT(self):
        if cluster_env == "development" or cluster_env == "test":
            return redis.Redis(host="host.docker.internal", port=6379)

        return redis.Redis(
            host=self.REDIS_ENDPOINT["Address"],
            port=self.REDIS_ENDPOINT["Port"],
            ssl=True,
            ssl_cert_reqs=None,
        )

if cluster_env == "development" or cluster_env == "test":
    aws_account_id = "000000000000"

config = Config(
    CLUSTER_ENV=cluster_env,
    TIMEOUT=timeout,
    IGNORE_TIMEOUT=ignore_timeout,
    AWS_ACCOUNT_ID=aws_account_id,
    AWS_DEFAULT_REGION=aws_region,
    BOTO_RESOURCE_KWARGS={"region_name": aws_region},
    CELL_SETS_BUCKET=f"cell-sets-{cluster_env}-{aws_account_id}",
    SOURCE_BUCKET=f"processed-matrix-{cluster_env}-{aws_account_id}",
    RESULTS_BUCKET=f"worker-results-{cluster_env}-{aws_account_id}",
    R_WORKER_URL="http://localhost:4000",
    # this works because in CI, `data/` is deployed under `worker/`
    # whereas in a container, it is mounted to `/data`. Either way, this ensures
    # that the appropriate path is selected, as both are two directories up
    LOCAL_DIR=os.path.join(os.pardir, os.pardir, "data"),
    RDS_PATH = "/data/processed.rds",
    TMP_RESULTS_PATH_GZ = "/data/rResult.gz"
)

config.API_URL = (
    f"http://api-{config.SANDBOX_ID}.api-{config.SANDBOX_ID}.svc.cluster.local:3000"
)

if cluster_env == "development" or cluster_env == "test":
    config.AWS_DEFAULT_REGION = "eu-west-1"
    config.BOTO_RESOURCE_KWARGS["aws_access_key_id"] = "my-key"
    config.BOTO_RESOURCE_KWARGS["aws_secret_access_key"] = "my-secret-key"
    config.BOTO_RESOURCE_KWARGS["region_name"] = "eu-west-1"
    config.API_URL = "http://host.docker.internal:3000"

    global_sdk_config.set_sdk_enabled(False)

if cluster_env == "development":
    config.BOTO_RESOURCE_KWARGS["endpoint_url"] = "http://host.docker.internal:4566"
    config.R_WORKER_URL = "http://r:4000"

if cluster_env != "test":
    core.patch_all()

if not config.AWS_ACCOUNT_ID:
    raise Exception("env.AWS_ACCOUNT_ID is not set")

if not config.AWS_DEFAULT_REGION:
    raise Exception("env.AWS_DEFAULT_REGION is not set")