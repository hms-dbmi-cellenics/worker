import datetime
import backoff
import requests
import os
from datetime import timezone
from logging import error, info

import aws_xray_sdk as xray
import boto3
from aws_xray_sdk.core import xray_recorder
from socket_io_emitter import Emitter

from ..config import config
from ..constants import DOWNLOAD_EXPERIMENT,LOAD_EXPERIMENT


class CountMatrix:
    def __init__(self):
        self.config = config
        self.local_path = os.path.join(
            self.config.LOCAL_DIR, self.config.EXPERIMENT_ID
        )
        self.s3 = boto3.client("s3", **self.config.BOTO_RESOURCE_KWARGS)

        self.last_fetch = None

    def get_objects(self):
        objects = self.s3.list_objects_v2(
            Bucket=self.config.SOURCE_BUCKET, Prefix=self.config.EXPERIMENT_ID
        )
        objects = objects.get("Contents")

        if not objects:
            return {}

        objects = {
            o["Key"]: o["LastModified"] for o in objects if o["Size"] > 0
        }

        return objects

    @xray_recorder.capture("CountMatrix.download_object")
    def download_object(self, key, last_modified):
        path = os.path.join(self.config.LOCAL_DIR, key)

        try:
            last_mod_local = datetime.datetime.fromtimestamp(
                os.path.getmtime(path), tz=timezone.utc
            )
        except FileNotFoundError:
            last_mod_local = None
        except Exception as e:
            error(e)
            last_mod_local = None

        if self.last_fetch and last_modified < self.last_fetch:
            info(
                f"Did not fetch as last modified (remote) of {last_modified}"
                f" was before last fetch time of {self.last_fetch}"
            )

            return False
        elif last_mod_local and last_modified < last_mod_local:
            info(
                f"Did not fetch as last modified (remote) of {last_modified}"
                f" was before last modified (local) of {last_mod_local}"
            )

            return False
        else:
            info(
                f"Fetching as last modified date of {last_modified}"
                f" is more recent than {self.last_fetch or 'Never'}"
            )

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        with open(path, "wb+") as f:
            io = Emitter({"client": config.REDIS_CLIENT})
            io.Emit(f'Heartbeat-{self.config.EXPERIMENT_ID}', {"type": "WorkResponse", "workingOn": DOWNLOAD_EXPERIMENT})
            info(f"Downloading {key} from S3...")
            self.s3.download_fileobj(
                Bucket=self.config.SOURCE_BUCKET,
                Key=key,
                Fileobj=f,
            )

            io.Emit(f'Heartbeat-{self.config.EXPERIMENT_ID}', {"type": "WorkResponse", "workingOn": LOAD_EXPERIMENT})

            self.last_fetch = last_modified
            f.seek(0)

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        return True

    @backoff.on_exception(
        backoff.constant, requests.exceptions.RequestException, interval=5
    )
    def check_if_received(self):
        info('Count matrices updated, checking if R worker is alive...')
        r = requests.get(
            f"{config.R_WORKER_URL}/health",
        )


    @xray_recorder.capture("CountMatrix.sync")
    def sync(self):

        # check if path existed before running this
        self.path_exists = os.path.exists(self.local_path)

        if not self.path_exists:
            info(f"Path {self.local_path} does not yet exist, creating it...")
            os.makedirs(self.local_path)

        objects = self.get_objects()

        info(f"Found {len(objects)} objects matching experiment.")
        synced = {
            key: self.download_object(key, last_modified)
            for key, last_modified in objects.items()
        }

        if True in synced.values():
            self.check_if_received()