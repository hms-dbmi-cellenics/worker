worker-python
=============

Overview
--------
The purpose of the worker is to carry out data analysis tasks. `worker-python` is the main entry point
for data analysis tasks. Each bioinformatics task to be performed is received from an AWS SQS queue.
Results are submitted to an SNS topic that is listened to by the `api` module.

Running the worker
------------------

### Start the worker

See the main README for instructions on how to run the workers in Docker.

### Process tasks
Tasks are automatically processed when they are received from the SQS queue specified.

For **local development**, make sure you have [InfraMock](https://github.com/biomage-ltd/inframock)
running alongside the [ui](https://github.com/biomage-ltd/ui) and [api](https://github.com/biomage-ltd/api)
projects. Refer to their respective documentations on how to run them locally. Once all of these are running,
tasks should automatically be submitted and processed when you perform actions on the `ui`. There is nothing else to do.

### Advanced: using the worker with live AWS queues

If you want to fetch and process data from **the live** `staging` or `production` clusters, you need to set
the `CLUSTER_ENV` environment variable to be either `staging` or `production`, respectively. The name of the SQS queue used
by the worker to subsribe to work is stored in the `WORKER_QUEUE` environment variable. You can set it to any active queue
that you know is currently being used, but it defaults to `development-queue.fifo`. This is created specifically and only
for local testing purposes.

In order to push work to the worker you can either:
* a) Submit work via the AWS console, using the SQS interface
* b) Submit work from the terminal using `aws-cli`

### Advanced: pushing custom work to the local worker

You can also push work to a locally running `worker` instance without using the `ui` or the `api` projects.

First, make sure you have InfraMock running. Then, you can use `aws-cli` to send a payload directly to the queue
the worker is listening to:

    aws --endpoint-url=http://localhost:4566 sqs send-message --queue-url http://localhost:4566/queue/development-queue.fifo --message-body "$(< payload.json)" --message-group-id "$(date +%s)"

This will push the payload in the file `payload.json` to the SQS queue the worker is listening to.


Development
-----------

Open the Visual Studio Code workspace:

    code python.code-workspace

You should be prompted to run the code in a container. If not, hit `Cmd+Shift+P` and search for
`Remote Containers: Reopen in Container`. Selecting this item will cause the Python worker to
reopen in a container for development.

### Tests
While having the workspace running in a container, open a terminal in VS code. Go to `src/` and run:

    gunzip -k ../../data/test/python.h5ad
    CLUSTER_ENV="test" python -m pytest --cov=.

### Task formatting
Task definitions are stored in the `api` project as an OpenAPI schema.
You can find this [here](https://github.com/biomage-ltd/api/blob/master/src/specs/api.yaml).

Download the schema and open it using Stoplight Studio. Looking into `WorkRequest` should give you
the schemas and parameters for all supported tasks.

### Sample tasks
Here are some examples:

* `GetEmbedding`:

    {
        "uuid": "509520fe-d329-437d-8752-b5868ad59425",
        "socketId": "Y1poEygzBfrDmIWpAAAA",
        "experimentId": "5928a56c7cbff9de78974ab50765ed20",
        "timeout": "2021-01-01T00:00:00Z",
        "body": {
            "name": "GetEmbedding",
            "type": "pca"
        }
    }

* `ListGenes`:

    {
        "uuid": "509520fe-d329-437d-8752-b5868ad59425",
        "socketId": "Y1poEygzBfrDmIWpAAAA",
        "experimentId": "5928a56c7cbff9de78974ab50765ed20",
        "timeout": "2021-01-01T00:00:00Z",
        "body": {
            "name": "ListGenes",
            "selectFields": ["highly_variable", "gene_names", "dispersions"],
            <!-- "geneNamesFilter": "%IN%", add this to filter results so only gene_names that contain IN in their names appear -->
            "orderBy": "dispersions",
            "orderDirection": "desc",
            "offset": 0,
            "limit": 20
        }
    }

* `GeneExpression`:

    {
        "uuid": "509520fe-d329-437d-8752-b5868ad59425",
        "socketId": "Y1poEygzBfrDmIWpAAAA",
        "experimentId": "5928a56c7cbff9de78974ab50765ed20",
        "timeout": "2099-12-31 00:00:00",
        "body": {
            "name": "GeneExpression",
            "cellSets": ["louvain-0", "louvain-1", "louvain-2", "louvain-3", "louvain-4", "louvain-5", "louvain-6"],
            "genes": ["TGFB1", "CST3"]
        }
    }

