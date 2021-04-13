import boto3
import datetime
import os
import hashlib
import requests
import backoff
from datetime import timezone
from config import get_config
import aws_xray_sdk as xray
from aws_xray_sdk.core import xray_recorder
import requests

config = get_config()


class CountMatrix:
    def __init__(self):
        self.config = get_config()
        self.local_path = os.path.join(self.config.LOCAL_DIR, self.config.EXPERIMENT_ID)
        self.s3 = boto3.client("s3", **self.config.BOTO_RESOURCE_KWARGS)

        self.last_fetch = None

    def get_objects(self):
        objects = self.s3.list_objects_v2(
            Bucket=self.config.SOURCE_BUCKET, Prefix=self.config.EXPERIMENT_ID
        )
        objects = objects.get("Contents")

        if not objects:
            return {}

        objects = {o["Key"]: o["LastModified"] for o in objects if o["Size"] > 0}

        return objects

    @xray_recorder.capture("CountMatrix.download_object")
    def download_object(self, key, last_modified):
        path = os.path.join(config.LOCAL_DIR, key)
        print(key)

        last_mod_local = None

        try:
            last_mod_local = datetime.datetime.fromtimestamp(
                os.path.getmtime(path), tz=timezone.utc
            )
        except Exception as e:
            print(e)
            last_mod_local = None

        if self.last_fetch and last_modified < self.last_fetch:
            print(
                datetime.datetime.utcnow(),
                "Did not fetch as last modified (remote) of",
                last_modified,
                "was before last fetch time of",
                self.last_fetch,
            )

            return False
        elif last_mod_local and last_modified < last_mod_local:
            print(
                datetime.datetime.utcnow(),
                "Did not fetch as last modified (remote) of",
                last_modified,
                "was before last modified (local) of",
                last_mod_local,
            )

            return False
        else:
            print(
                datetime.datetime.utcnow(),
                "Fetching as last modified date of",
                last_modified,
                "is more recent than",
                self.last_fetch or "Never",
            )

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        with open(path, "wb+") as f:
            self.s3.download_fileobj(
                Bucket=self.config.SOURCE_BUCKET,
                Key=key,
                Fileobj=f,
            )

            self.last_fetch = last_modified
            f.seek(0)

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        return True

    @xray_recorder.capture("CountMatrix.sync")
    def sync(self):
        # check if path existed before running this
        self.path_exists = os.path.exists(self.local_path)

        if not self.path_exists:
            print(
                datetime.datetime.utcnow(),
                "Path",
                self.local_path,
                "does not yet exist, creating it...",
            )
            os.makedirs(self.local_path)

        objects = self.get_objects()

        print(
            datetime.datetime.utcnow(),
            "Found",
            len(objects),
            "objects matching experiment.",
        )

        synced = {
            key: self.download_object(key, last_modified)
            for key, last_modified in objects.items()
        }
