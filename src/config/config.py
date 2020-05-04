import os


class BaseConfig(object):
    ENVIRONMENT = "base"
    DYNAMO_TABLE = "experiments-staging"
    QUEUE_NAME = os.getenv("WORK_QUEUE", default="test-queue")
    RESULTS_BUCKET = "worker-results-staging"
    TIMEOUT = int(os.getenv("WORK_TIMEOUT", default="1200"))


class TestConfig:
    ENVIRONMENT = "testing"


class DevelopmentConfig(BaseConfig):
    ENVIRONMENT = "development"


class StagingConfig(BaseConfig):
    ENVIRONMENT = "staging"


class ProductionConfig(StagingConfig):
    ENVIRONMENT = "production"
