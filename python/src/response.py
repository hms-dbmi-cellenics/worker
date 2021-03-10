import boto3
import json
import datetime
from functools import reduce
from config import get_config
import uuid

config = get_config()


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

    def _get_response_msg(self, s3_keys=None):
        if s3_keys:
            result_objs = []

            for (res, s3_key) in zip(self.results, s3_keys):
                s3_path = "/".join([self.s3_bucket, s3_key])

                obj = res.get_result_object(resp_format=True, s3_path=s3_path)
                result_objs.append(obj)
        else:
            result_objs = [
                res.get_result_object(resp_format=True) for res in self.results
            ]

        return {
            "request": self.request,
            "results": list(result_objs),
            "response": {"cacheable": self.cacheable, "error": self.error},
        }

    def _upload(self, result):
        client = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        key = "{}/{}".format(self.request["uuid"], str(uuid.uuid4()))
        body = result.get_result_object()["body"]

        client.put_object(Key=key, Bucket=self.s3_bucket, Body=body)
        return key

    def _send_notification(self, mssg):
        sns = boto3.client("sns", **config.BOTO_RESOURCE_KWARGS)

        msg_to_send = json.dumps({"default": json.dumps(mssg)})

        r = sns.publish(
            TargetArn="arn:aws:sns:{}:{}:{}".format(
                config.AWS_REGION, config.AWS_ACCOUNT_ID, config.SNS_TOPIC
            ),
            Message=msg_to_send,
            MessageAttributes={
                'type': {
                    'DataType': 'String',
                    'StringValue': 'WorkResponse'
                }
            },
            MessageStructure="json",
        )

        print(datetime.datetime.utcnow(), "Message successfully sent to sns", r)

        return msg_to_send

    def publish(self):
        # Get total length of all result objects:
        message_length = reduce(
            lambda acc, curr: acc + curr.get_result_length(), self.results, 0,
        ) + len(json.dumps(self.request))

        # If we are over 80% of the limit (256 KB, 262144 bytes), upload to S3.
        # Otherwise, we can send the entire payload through the SNS topic.
        MAX_SNS_MESSAGE_LEN = 262144
        upload_to_s3 = message_length >= 0.8 * MAX_SNS_MESSAGE_LEN

        if upload_to_s3:
            s3_keys = [self._upload(res) for res in self.results]
            response_msg = self._get_response_msg(s3_keys=s3_keys)
        else:
            response_msg = self._get_response_msg()

        print(response_msg)

        return self._send_notification(response_msg)
