import os
from unittest.mock import mock_open, patch

from worker.config import config


class TestExperimentIDFetch:
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data='key="value"\ntest="test"\nexperimentId="mockExperiment"\nsandboxId="mockSandbox"\n',
    )
    def test_config_reads_labels_from_file(self, mock_file):
        assert config.EXPERIMENT_ID == "mockExperiment"
        assert config.SANDBOX_ID == "mockSandbox"

    def test_config_reads_labels_from_env_if_file_not_found(self):
        os.environ["EXPERIMENT_ID"] = "fake_experiment"
        os.environ["SANDBOX_ID"] = "fake_sandbox"

        assert config.EXPERIMENT_ID == os.getenv("EXPERIMENT_ID")
        assert config.SANDBOX_ID == os.getenv("SANDBOX_ID")

    def test_config_returns_none_when_no_data_found(self):

        try:
            del os.environ["EXPERIMENT_ID"]
            del os.environ["SANDBOX_ID"]
        except KeyError:
            pass

        assert config.EXPERIMENT_ID is None
        assert config.SANDBOX_ID is None
