import os

from .config import (
    DevelopmentConfig,
    StagingConfig,
    ProductionConfig,
)

config = {
    "default": DevelopmentConfig,
    "staging": StagingConfig,
    "production": ProductionConfig,
}


def get_config():
    env_name = os.getenv("GITLAB_ENVIRONMENT_NAME", "default")

    config_class = config.get(env_name, config["default"])
    return config_class()
