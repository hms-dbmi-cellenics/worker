import os

from .config import (
    TestConfig,
    DevelopmentConfig,
    StagingConfig,
    ProductionConfig,
)

config = {
    "default": DevelopmentConfig,
    "dev": DevelopmentConfig,
    "testing": TestConfig,
    "staging": StagingConfig,
    "production": ProductionConfig,
}


def get_config():
    env_name = os.getenv("GITLAB_ENVIRONMENT_NAME", "default")
    return config[env_name]
