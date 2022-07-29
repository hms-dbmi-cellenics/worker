import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..result import Result
from ..tasks import Task
from worker.config import config

import PIL
from PIL import Image
import numpy as np

import boto3
s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)

class GetImgPlot(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.task_etag = msg["ETag"]
        self.experiment_id=msg["experimentId"]
        self.name = 'GetImgPlot'

    def _format_result(self, result):
        # Return a list of formatted results.
        return Result(result, upload=False)

    def _format_request(self):
        request = {
            "plotType": self.task_def["type"],
            "genes": self.task_def["features"],
            "etag": self.task_etag,
            "plotSubType": self.task_def["plotSubType"],
        }
        return request

    def _generate_img(self, pixels):
        rmat = pixels[0]
        gmat = pixels[1]
        bmat = pixels[2]
        # if len(pixels) == 4:
        #     amat = pixels[3]
        # else:
                
        formatted_pixels = []
        for rrow, grow, brow in zip(rmat, gmat, bmat):
            formatted_row = [(r, g, b) for r, g, b in zip(rrow, grow, brow)]
            formatted_pixels.append(formatted_row)

        ready = np.array(formatted_pixels, dtype=np.uint8)
        img = Image.fromarray(ready)
        img.save('img.png')
        s3.upload_file(
            Filename = './img.png',
            Bucket = config.RESULTS_BUCKET,
            Key = self.task_etag
        )
        print('EXP ID ', self.experiment_id, 'ETAG ', self.task_etag, 'name ', self.name)
        s3.put_object_tagging(
            Key=self.task_etag,
            Bucket=config.RESULTS_BUCKET,
            Tagging={
                "TagSet": [
                    {"Key": "experimentId", "Value": self.experiment_id},
                    {"Key": "requestType", "Value": self.name},
                ]
            },
        )

    @xray_recorder.capture("GetImgPlot.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()
        response = requests.post(
            f"{config.R_WORKER_URL}/v0/getImgPlot",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)
        raw_data = result.get("data")
        obj = json.loads(raw_data)
        
        data = obj['data']
        
        self._generate_img(data)
        
        return self._format_result({})
