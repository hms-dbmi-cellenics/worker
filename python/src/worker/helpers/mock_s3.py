import json
from operator import itemgetter

cell_set_responses = {
    "one_set": [
        {
            "name": "my amazing cluster",
            "key": "cluster1",
            "cellIds": [4, 5],
        },
    ],
    "two_sets": [
        {
            "name": "my amazing cluster",
            "key": "cluster1",
            "cellIds": [4, 5],
        },
        {
            "name": "my other amazing cluster",
            "key": "cluster2",
            "cellIds": [0, 1, 2, 3],
        },
    ],
    "two_sets_intersected": [
        {
            "name": "intersecting set",
            "key": "cluster1",
            "cellIds": [1, 2, 3],
        },
        {
            "name": "other intersecting set",
            "key": "cluster2",
            "cellIds": [3, 4, 5],
        },
    ],
    "three_sets": [
        {
            "name": "one set",
            "key": "cluster1",
            "cellIds": [4, 5],
        },
        {
            "name": "other set",
            "key": "cluster2",
            "cellIds": [0, 1, 2, 3],
        },
        {
            "name": "basis set",
            "key": "basisCluster",
            "cellIds": [0, 1, 5],
        },
    ],
    "hierarchichal_sets": [
        {
            "name": "hierarchy 1",
            "key": "set_hierarchy_1",
            "cellIds": [],
            "children": [
                {
                    "name": "one set",
                    "key": "cluster1",
                    "cellIds": [4],
                },
                {
                    "name": "another set",
                    "key": "cluster2",
                    "cellIds": [5, 6],
                },
            ],
        },
        {
            "name": "hierarchy 2",
            "key": "set_hierarchy_2",
            "cellIds": [],
            "children": [
                {
                    "name": "set",
                    "key": "cluster3",
                    "cellIds": [0, 1, 2, 3],
                },
                {
                    "name": "set1",
                    "key": "cluster4",
                    "cellIds": [1, 10, 11, 12, 13, 14],
                },
            ],
        },
    ],
}


class MockS3Class:
    response = cell_set_responses["one_set"]

    def setResponse(response_key):
        MockS3Class.response = cell_set_responses[response_key]

    def download_fileobj(*args, **kwargs):
        Bucket, Key, Fileobj = itemgetter("Bucket", "Key", "Fileobj")(kwargs)

        if not Bucket or not Key or not Fileobj:
            raise Exception("Parameters not received")

        Fileobj.write(str.encode(f'{{"cellSets": {json.dumps(MockS3Class.response)}}}'))

        return
