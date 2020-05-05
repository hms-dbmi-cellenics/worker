import boto3
import json
import datetime

from config import get_config

config = get_config()


class Result:
    def __init__(self, work_def, result):
        self.uuid = work_def["uuid"]
        self.socket_id = work_def["socketId"]
        self.result = result.tolist()
        self.s3_bucket = config.RESULTS_BUCKET
        self.s3_key = self.uuid

    def _get_response_schema(self):
        return {
            "uuid": self.uuid,
            "socketId": self.socket_id,
            "timeout": "2021-09-09T12:00:00Z",
            "results": [
                {
                    "content-type": "application/json",
                    "type": "s3-path",
                    "body": "/".join([self.s3_bucket, self.s3_key]),
                }
            ],
        }

    def _upload(self):
        client = boto3.client("s3")
        client.put_object(
            Key=self.s3_key, Bucket=self.s3_bucket, Body=json.dumps(self.result)
        )
        print(datetime.datetime.now(), "Result was successfully uploaded to s3.")

    def _send_notification(self):
        mssg = self._get_response_schema()
        sns = boto3.client("sns")
        sns.publish(
            TargetArn="arn:aws:sns:eu-west-2:242905224710:work-results-staging",
            Message=json.dumps({"default": json.dumps(mssg)}),
            MessageStructure="json",
        )
        print(
            datetime.datetime.now(), "Message {} successfully sent to sns".format(mssg)
        )

    def publish(self):
        self._upload()
        self._send_notification()
        print(datetime.datetime.now(), "Result was published.")
