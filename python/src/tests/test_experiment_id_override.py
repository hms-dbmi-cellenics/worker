import os
from unittest.mock import mock_open, patch

from worker.config import config


class TestExperimentIDFetch:
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data='key="value"\nsandboxId="mockSandbox"\n',
    )
    def test_config_reads_labels_from_file(self, mock_file):
        assert config.SANDBOX_ID == "mockSandbox"

    def test_config_reads_labels_from_env_if_file_not_found(self):
        os.environ["SANDBOX_ID"] = "fake_sandbox"

        assert config.SANDBOX_ID == os.getenv("SANDBOX_ID")

    def test_config_returns_none_when_no_data_found(self):
        try:
            del os.environ["SANDBOX_ID"]
        except KeyError:
            pass

        assert config.SANDBOX_ID is None
