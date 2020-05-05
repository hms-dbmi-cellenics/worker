import boto3
import json
import datetime

from config import get_config

config = get_config()


class Result:
    def __init__(self, work_def, result):
        self.uuid = work_def["uuid"]
        self.socket_id = work_def["socketId"]
        self.result = result
        self.s3_bucket = config.RESULTS_BUCKET
        self.s3_key = self.uuid

    def _get_response_mssg(self, is_s3):
        if is_s3:
            result_type = "s3-path"
            result_body = "/".join([self.s3_bucket, self.s3_key])
        else:
            result_type = "inline"
            result_body = str(self.result)

        return {
            "uuid": self.uuid,
            "socketId": self.socket_id,
            "timeout": "2021-09-09T12:00:00Z",
            "results": [
                {
                    "content-encoding": "utf-8",
                    "content-type": "application/json",
                    "type": result_type,
                    "body": result_body,
                }
            ],
        }

    def _upload(self):
        client = boto3.client("s3")
        client.put_object(
            Key=self.s3_key, Bucket=self.s3_bucket, Body=json.dumps(self.result)
        )
        print(datetime.datetime.now(), "Result was successfully uploaded to s3.")

    def _send_notification(self, mssg):
        sns = boto3.client("sns")
        sns.publish(
            TargetArn="arn:aws:sns:{}:{}:{}".format(
                config.AWS_REGION, config.AWS_ACCOUNT_ID, config.SNS_TOPIC
            ),
            Message=json.dumps({"default": json.dumps(mssg)}),
            MessageStructure="json",
        )
        print(datetime.datetime.now(), "Message successfully sent to sns")

    def publish(self):
        # if the result is too big for SNS (bigger than 256KB), put it in S3 first. Note: I leave 4000 chars for the message schema
        chars_numb = len(str(self.result))
        put_in_s3 = chars_numb / 1000 >= 252

        if put_in_s3:
            self._upload()

        mssg = self._get_response_mssg(put_in_s3)
        self._send_notification(mssg)
        print(datetime.datetime.now(), "Result was published.")
