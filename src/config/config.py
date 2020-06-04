import os

"""
NOTE: We only have two deployments at the moment, `staging` and `production`.

For this reason, all environments that are not production are going to
use the staging deployments.

This will change once we have dedicated development/staging/review/etc.
environments.
"""


class BaseConfig(object):
    ENVIRONMENT = "development"

    # When we are running e.g. review deployments for MRs, we want
    # to use the `staging` deployment for all other services, because
    # we do not have a dedicated testing/dev environment right now.
    # This will need to change as those get deployed eventually.
    CLUSTER_ENVIRONMENT = "staging"

    QUEUE_NAME = os.getenv("WORK_QUEUE", default="test-queue")
    TIMEOUT = int(os.getenv("WORK_TIMEOUT", default="1200"))
    AWS_ACCOUNT_ID = "242905224710"
    AWS_REGION = "eu-west-2"

    def get_dynamo_table(self):
        return f"experiments-{self.CLUSTER_ENVIRONMENT}"

    def get_results_bucket(self):
        return f"worker-results-{self.CLUSTER_ENVIRONMENT}"

    def get_sns_topic(self):
        return f"work-results-{self.CLUSTER_ENVIRONMENT}"


class DevelopmentConfig(BaseConfig):
    ENVIRONMENT = "development"
    CLUSTER_ENVIRONMENT = "staging"


class StagingConfig(BaseConfig):
    ENVIRONMENT = "staging"


class ProductionConfig(StagingConfig):
    ENVIRONMENT = "production"
    CLUSTER_ENVIRONMENT = "production"
