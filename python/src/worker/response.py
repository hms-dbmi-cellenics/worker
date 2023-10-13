import gzip
import io
import ujson
from logging import info
import base64
import sys

import aws_xray_sdk as xray
import boto3
from aws_xray_sdk.core import xray_recorder
from socket_io_emitter import Emitter
from io import BytesIO

from .config import config
from worker_status_codes import (
    COMPRESSING_TASK_DATA,
    UPLOADING_TASK_DATA,
    FINISHED_TASK,
)
from worker.helpers.send_status_updates import send_status_update

from datetime import datetime

class Response:
    def __init__(self, request, result):
        self.request = request
        self.result = result

        self.error = result.error
        self.cacheable = (not result.error) and result.cacheable

        self.s3_bucket = config.RESULTS_BUCKET

    #' Returns the compressed work result to be sent
    #'
    #' @return gzipped_body to upload to s3 and the compressed bytes 
    #' object to send over redis if the work result is small enough
    def _construct_data_for_upload(self):
        info("Starting compression before upload to s3")
        io = Emitter({"client": config.REDIS_CLIENT})
        send_status_update(
            io, self.request["experimentId"], COMPRESSING_TASK_DATA, self.request
        )

        gzipped_body = BytesIO()
        with gzip.open(gzipped_body, "wt", encoding="utf-8") as zipfile:
            if isinstance(self.result.data, str):
                info("Compressing string work result")
                zipfile.write(self.result.data)
            else:
                info('Encoding and compressing json work result')
                ujson.dump(self.result.data, zipfile)

        gz_body_bytes = None
        # If size is less than 250 kb, then send it over notification too
        # Needs to be done here because upload_fileobj closes the file:
        # https://github.com/boto/boto3/issues/929
        kb = 1000
        body_size = sys.getsizeof(gzipped_body)
        info(f"Body size is {body_size}")
        if (body_size <= 250 * kb):
            info("Data is smaller than 250 kb, sending over socket")
            gzipped_body.seek(0)
            gz_body_bytes = gzipped_body.read()

        gzipped_body.seek(0)

        info("Compression finished")
        return gzipped_body, gz_body_bytes

    def _construct_response_msg(self, socket_data = None):
        if socket_data:
            return base64.b64encode(socket_data)

        message = {
            "request": self.request,
            "response": {"cacheable": self.cacheable, "error": self.error, "signedUrl": self.request["signedUrl"]},
            "type": "WorkResponse",
        }

        if self.error:
            message["response"]["errorCode"] = self.result.data["error_code"]
            message["response"]["userMessage"] = self.result.data["user_message"]

        return message

    @xray_recorder.capture("Response._upload")
    def _upload(self, response_data, type):
        io = Emitter({"client": config.REDIS_CLIENT})
        send_status_update(
            io, self.request["experimentId"], UPLOADING_TASK_DATA, self.request
        )

        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        ETag = self.request["ETag"]

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        if type == "path":
            with open(response_data, "rb") as file:
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

    #' Send a notification that a work response finished
    #'
    #' @param socket_data Optional. The work result, if not None, it is sent instead
    #'  of the default json response msg so it reaches the client faster
    #'
    #' @export
    def _send_notification(self, socket_data=None):
        io = Emitter({"client": config.REDIS_CLIENT})
        if self.request.get("broadcast"):
            io.Emit(
                f'ExperimentUpdates-{self.request["experimentId"]}',
                self._construct_response_msg(),
            )

            info(f"Broadcast results to users viewing experiment {self.request['experimentId']}.")

        send_status_update(
            io, self.request["experimentId"], FINISHED_TASK, self.request
        )

        io.Emit(f'WorkResponse-{self.request["ETag"]}', self._construct_response_msg(socket_data))

        info(f"Notified users waiting for request with ETag {self.request['ETag']}.")

    @xray_recorder.capture("Response.publish")
    def publish(self):
        info(f"Request {self.request['ETag']} processed, response:")

        socket_data = None

        if not self.error and self.cacheable:
            info("Uploading response to S3")
            if self.result.data == config.RDS_PATH or self.result.data == config.INTERNAL_RESULTS_PATH:
                self._upload(self.result.data, "path")
            else:
                s3_data, socket_data = self._construct_data_for_upload()
                self._upload(s3_data, "obj")

        info("Sending socket.io message to clients subscribed to work response")
        return self._send_notification(socket_data)
