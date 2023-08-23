import gzip
import io
import json
from logging import info

import aws_xray_sdk as xray
import boto3
from aws_xray_sdk.core import xray_recorder
from socket_io_emitter import Emitter

from .config import config
from .constants import COMPRESSING_TASK_DATA, UPLOADING_TASK_DATA, FINISHED_TASK

class Response:
    def __init__(self, request, result):
        self.request = request
        self.result = result

        self.error = result.error
        self.cacheable = (not result.error) and result.cacheable

        self.s3_bucket = config.RESULTS_BUCKET

    def _construct_data_for_upload(self):
        info("Starting compression before upload to s3")
        gzipped_body = io.BytesIO()
        with gzip.open(gzipped_body, "wt", encoding="utf-8") as zipfile:
            if (isinstance(self.result.data, str)):
                info('Compressing string work result')
                zipfile.write(self.result.data)
            else:
                info('Encoding and compressing json work result')
                json.dump(self.result.data, zipfile)


        gzipped_body.seek(0)

        info("Compression finished")
        return gzipped_body

    def _construct_response_msg(self):
        message = {
            "request": self.request,
            "response": {"cacheable": self.cacheable, "error": self.error},
            "type": "WorkResponse",
        }

        if self.error:
            message["response"]["errorCode"] = self.result.data["error_code"]
            message["response"]["userMessage"] = self.result.data["user_message"]

        return message

    @xray_recorder.capture("Response._upload")
    def _upload(self, response_data, type):
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        ETag = self.request["ETag"]

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        if (type == "path"):
            with open(response_data, 'rb') as file:
                client.upload_fileobj(file, self.s3_bucket, ETag)
        else:
            client.upload_fileobj(response_data, self.s3_bucket, ETag)


        client.put_object_tagging(
            Key=ETag,
            Bucket=self.s3_bucket,
            Tagging={
                "TagSet": [
                    {"Key": "experimentId", "Value": self.request["experimentId"]},
                    {"Key": "requestType", "Value": self.request["body"]["name"]},
                ]
            },
        )

        info(f"Response was uploaded in bucket {self.s3_bucket} at key {ETag}.")

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        return ETag

    def _send_notification(self):
        io = Emitter({"client": config.REDIS_CLIENT})
        if self.request.get("broadcast"):
            io.Emit(
                f'ExperimentUpdates-{self.request["experimentId"]}',
                self._construct_response_msg(),
            )

            info(
                f"Broadcast results to users viewing experiment {self.request['experimentId']}."
            )

        io.Emit(f'Heartbeat-{self.request["experimentId"]}', {"type": "WorkResponse", "workingOn": FINISHED_TASK, "request": self.request})
        io.Emit(f'WorkResponse-{self.request["ETag"]}', self._construct_response_msg())

        info(f"Notified users waiting for request with ETag {self.request['ETag']}.")

    @xray_recorder.capture("Response.publish")
    def publish(self):
        info(f"Request {self.request['ETag']} processed, response:")

        if not self.error and self.cacheable:
            io = Emitter({"client": config.REDIS_CLIENT})

            io.Emit(f'Heartbeat-{self.request["experimentId"]}',
             {"type": "WorkResponse", "workingOn": COMPRESSING_TASK_DATA, "request": self.request})

            response_data = self._construct_data_for_upload()

            io.Emit(f'Heartbeat-{self.request["experimentId"]}',
             {"type": "WorkResponse", "workingOn": UPLOADING_TASK_DATA, "request": self.request})
            info("Uploading response to S3")
            if (self.result.data == config.RDS_PATH):
                response_data = self.result.data
                self._upload(response_data, "path")
            else:
                response_data = self._construct_data_for_upload()
                self._upload(response_data, "obj")

        info("Sending socket.io message to clients subscribed to work response")
        return self._send_notification()
