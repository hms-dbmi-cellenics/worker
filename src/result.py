import boto3
import json
import os


class Result:
    def __init__(self, work_def, result):
        self.uuid = work_def.get("uuid", "1234")
        self.socket_id = work_def.get("socketId", "567")
        self.result = {"result": result.tolist()}
        self.s3_bucket = os.getenv("RESULTS_BUCKET", default="worker-results-staging")
        self.s3_key = self.uuid

    def _get_response_schema(self):
        return {
            "uuid": self.uuid,
            "socketId": self.socket_id,
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
        print("Result was successfully uploaded to s3.")

    def _send_notification(self):
        mssg = self._get_response_schema()
        sns = boto3.client("sns")
        sns.publish(
            TargetArn="arn:aws:sns:eu-west-2:242905224710:work-results-staging",
            Message=json.dumps({"default": json.dumps(mssg)}),
            MessageStructure="json",
        )
        print("Message {} successfully sent to sns".format(mssg))

    def publish(self):
        self._upload()
        self._send_notification()
        print("Result was published.")
