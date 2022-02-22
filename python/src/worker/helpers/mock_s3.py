import json
import os
from operator import itemgetter

cell_set_responses = {}
with open(os.path.join("tests/data", "cell_sets.json")) as f:
    cell_set_responses = json.load(f)


class MockS3Class:
    response = cell_set_responses["one_set"]

    def setResponse(response_key):
        MockS3Class.response = cell_set_responses[response_key]

    def download_fileobj(*args, **kwargs):
        Bucket, Key, Fileobj = itemgetter("Bucket", "Key", "Fileobj")(kwargs)

        if not Fileobj:
            raise Exception("Parameters not received")

        Fileobj.write(str.encode(f'{{"cellSets": {json.dumps(MockS3Class.response)}}}'))

        return
