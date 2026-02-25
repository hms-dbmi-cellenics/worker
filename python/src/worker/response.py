import gzip
import orjson
import ujson
import zstandard as zstd
from logging import info
import base64
import sys
import os

import time
from datetime import datetime


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


def _format_bytes(num_bytes):
    """Convert bytes to human-readable format (B, KB, MB, GB)"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if num_bytes < 1000:
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1000
    return f"{num_bytes:.1f} TB"

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
    def _construct_data_for_upload_OLD(self):
        """Old version for comparison"""
        overall_start = time.time()
        
        info("Starting compression before upload to s3 (OLD METHOD)")
        
        serialize_start = time.time()
        gzipped_body = BytesIO()
        with gzip.open(gzipped_body, "wt", encoding="utf-8") as zipfile:
            if isinstance(self.result.data, str):
                zipfile.write(self.result.data)
            else:
                ujson.dump(self.result.data, zipfile)
        
        serialize_compress_time = time.time() - serialize_start
        compressed_size = gzipped_body.tell()
        info(f"⏱️  OLD: Serialize + Compress: {serialize_compress_time:.3f}s ({_format_bytes(compressed_size)})")

        gz_body_bytes = None
        kb = 1000
        body_size = sys.getsizeof(gzipped_body)
        
        if body_size <= 250 * kb:
            gzipped_body.seek(0)
            gz_body_bytes = gzipped_body.read()

        gzipped_body.seek(0)
        
        overall_time = time.time() - overall_start
        info(f"⏱️  OLD: Total: {overall_time:.3f}s")
        
        return gzipped_body, gz_body_bytes
    
    def _construct_data_for_upload(self):
        overall_start = datetime.now()
        
        info("Starting compression before upload to s3")
        io = Emitter({"client": config.REDIS_CLIENT})
        send_status_update(
            io, self.request["experimentId"], COMPRESSING_TASK_DATA, self.request
        )

        # Serialization timing
        serialize_start = datetime.now()
        if isinstance(self.result.data, str):
            info("Compressing string work result")
            data_bytes = self.result.data.encode('utf-8')
        else:
            info('Encoding and compressing json work result')
            data_bytes = orjson.dumps(self.result.data)
        
        serialize_time = (datetime.now() - serialize_start).total_seconds()
        uncompressed_size = len(data_bytes)
        info(f"⏱️  Serialization: {serialize_time:.3f}s ({_format_bytes(uncompressed_size)})")
        
        # Compression timing (zstd)
        compress_start = datetime.now()
        cctx = zstd.ZstdCompressor(level=3)
        compressed = cctx.compress(data_bytes)
        
        compress_time = (datetime.now() - compress_start).total_seconds()
        compressed_size = len(compressed)
        compression_ratio = (1 - compressed_size / uncompressed_size) * 100 if uncompressed_size > 0 else 0
        info(f"⏱️  Compression (zstd level 3): {compress_time:.3f}s ({_format_bytes(compressed_size)}, {compression_ratio:.1f}% reduction)")
        
        compressed_body = BytesIO(compressed)
        
        # Check if small enough to send over socket
        gz_body_bytes = None
        kb = 1000
        body_size = len(compressed)
        
        if body_size <= 250 * kb:
            info("Data is smaller than 250 kb, sending over socket")
            socket_start = datetime.now()
            gz_body_bytes = compressed
            socket_time = (datetime.now() - socket_start).total_seconds()
            info(f"⏱️  Socket read: {socket_time:.3f}s")

        compressed_body.seek(0)
        
        overall_time = (datetime.now() - overall_start).total_seconds()
        info(f"⏱️  Total compression (zstd): {overall_time:.3f}s")
        info("Compression finished")
        
        return compressed_body, gz_body_bytes

    def _construct_data_for_upload_NEW(self):
        overall_start = time.time()
        
        info("Starting compression before upload to s3")
        io = Emitter({"client": config.REDIS_CLIENT})
        send_status_update(
            io, self.request["experimentId"], COMPRESSING_TASK_DATA, self.request
        )

        # Serialization timing
        serialize_start = time.time()
        if isinstance(self.result.data, str):
            info("Compressing string work result")
            data_bytes = self.result.data.encode('utf-8')
        else:
            info('Encoding and compressing json work result')
            data_bytes = orjson.dumps(self.result.data)
        
        serialize_time = time.time() - serialize_start
        uncompressed_size = len(data_bytes)
        info(f"⏱️  Serialization: {serialize_time:.3f}s ({_format_bytes(uncompressed_size)})")
        
        # Compression timing
        compress_start = time.time()
        gzipped_body = BytesIO()
        with gzip.open(gzipped_body, "wb", compresslevel=6) as zipfile:
            zipfile.write(data_bytes)
        
        compress_time = time.time() - compress_start
        compressed_size = gzipped_body.tell()
        compression_ratio = (1 - compressed_size / uncompressed_size) * 100 if uncompressed_size > 0 else 0
        info(f"⏱️  Compression: {compress_time:.3f}s ({_format_bytes(compressed_size)}, {compression_ratio:.1f}% reduction)")

        # Check if small enough to send over socket
        gz_body_bytes = None
        kb = 1000
        body_size = gzipped_body.tell()
        
        if body_size <= 250 * kb:
            info("Data is smaller than 250 kb, sending over socket")
            socket_start = time.time()
            gzipped_body.seek(0)
            gz_body_bytes = gzipped_body.read()
            socket_time = time.time() - socket_start
            info(f"⏱️  Socket read: {socket_time:.3f}s")

        gzipped_body.seek(0)
        
        overall_time = time.time() - overall_start
        info(f"⏱️  Total compression: {overall_time:.3f}s")
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
        if self.request["requestProps"].get("broadcast"):
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
            if self.result.data == config.RDS_PATH or self.result.data == config.TMP_RESULTS_PATH_GZ:
                self._upload(self.result.data, "path")
            else:
                s3_data, socket_data = self._construct_data_for_upload()
                self._upload(s3_data, "obj")

        info("Sending socket.io message to clients subscribed to work response")
        self._send_notification(socket_data)

        # Remove the temporary file to transfer data between R and python
        # to free up memory
        if self.result.data == config.TMP_RESULTS_PATH_GZ:
            info("Cleaning up temp files generated by work result")
            os.remove(config.TMP_RESULTS_PATH_GZ)
