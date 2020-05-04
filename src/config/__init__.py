import os

from .config import (
    DevelopmentConfig,
    StagingConfig,
    ProductionConfig,
)

config = {
    "default": DevelopmentConfig,
    "dev": DevelopmentConfig,
    "staging": StagingConfig,
    "production": ProductionConfig,
}


def get_config():
    env_name = os.getenv("GITLAB_ENVIRONMENT_NAME", "default")
    return config[env_name]
