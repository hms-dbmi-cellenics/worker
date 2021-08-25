import json
import uuid
from functools import reduce
from logging import info

import aws_xray_sdk as xray
import boto3
from aws_xray_sdk.core import xray_recorder
from socket_io_emitter import Emitter

from .config import config


class Response:
    def __init__(self, request, results):
        self.request = request
        self.results = results

        self.error = False
        for result in self.results:
            self.error = self.error or result.error

        if self.error:
            self.cacheable = False
        else:
            self.cacheable = True
            for result in self.results:
                self.cacheable = self.cacheable and result.cacheable

        self.s3_bucket = config.RESULTS_BUCKET

    def _construct_response_msg(self, brief=False):
        message = {
            "request": self.request,
            "response": {"cacheable": self.cacheable, "error": self.error},
        }

        if not brief:
            result_objs = [
                res.get_result_object(resp_format=True) for res in self.results
            ]
            message["results"] = list(result_objs)
        
        return message

    @xray_recorder.capture("Response._upload")
    def _upload(self, response_msg):
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        ETag = self.request["ETag"]
        body = json.dumps(response_msg)

        # Disabled X-Ray to fix a botocore bug where the context
        # does not propagate to S3 requests. see:
        # https://github.com/open-telemetry/opentelemetry-python-contrib/issues/298
        was_enabled = xray.global_sdk_config.sdk_enabled()
        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(False)

        client.put_object(Key=ETag, Bucket=self.s3_bucket, Body=body)

        client.put_object_tagging(
            Key=ETag,
            Bucket=self.s3_bucket,
            Tagging={
                'TagSet': [
                    {
                        'Key': 'experimentId',
                        'Value': self.request['experimentId']
                    },
                    {
                        'Key': 'requestType',
                        'Value': self.request['body']['name']
                    },

                    # TODO: this needs to be removed and a proper
                    # ACL system implemented before merging.
                    {
                        'Key': 'public',
                        'Value': 'true',
                    },
                ]
            }
        )

        info(f"Repsonse was uploaded in bucket {self.s3_bucket} at key {ETag}.")

        if was_enabled:
            xray.global_sdk_config.set_sdk_enabled(True)

        return ETag

    def _send_notification(self, response_msg):
        io = Emitter({'host': config.REDIS_HOST, 'port':6379})

        if self.request["socketId"] == "broadcast":
            print(f'{self.request["experimentId"]}-{self.request["body"]["name"]}');

            io.Emit(
                f'{self.request["experimentId"]}-{self.request["body"]["name"]}',
                self._construct_response_msg(brief=True)
            )
        else:
            io.Emit(
                f'WorkResponse-{self.request["ETag"]}',
                self._construct_response_msg(brief=True)
            )

    
        info(f"Notified users waiting for request with ETag {self.request['ETag']}.")

    @xray_recorder.capture("Response.publish")
    def publish(self):
        response_msg = self._construct_response_msg()

        info(f"Request {self.request['ETag']} processed, response:")

        if not response_msg["response"]["error"] and response_msg["response"]["cacheable"]:
            info("Uploading response to S3")
            self._upload(response_msg)

        info("Sending socket.io message to clients subscribed to work response")
        return self._send_notification(response_msg)
