# Overview
The purpose of the worker is to carry out data analysis tasks, for example computing an embedding. Data analysis tasks get received from a SQS queue and results will get submitted to a SNS topic the API is subscribed to.


# Build and Deploy step
The deployment step for this repository consits of deploying all resources needed for work to be carried out. Note that the deployment step in this repository does not start the work. Work gets started dynamically upon a user request by the [API](https://github.com/biomage-ltd/api).

CI deployment needs four variables that need to be created on the CI/CD side:

* `GITLAB_DEPLOY_USER` and `GITLAB_DEPLOY_PASSWORD`. These need to specify a valid GitLab deploy key that has read rights to the Docker registry. This secret is used by dynamically generated jobs to automatically be able to pull the worker image.
* `K8S_SECRET_AWS_ACCESS_KEY_ID`, `K8S_SERET_AWS_SECRET_ACCESS_KEY`, `K8S_SECRET_AWS_DEFAULT_REGION`, which specify the IAM account and default region to be used with this worker.

# Development

## Run locally

### 1. Set up environment

See the main README on how to run the worker containers in docker-compose.

### 2. Submit a task
* If you want to fetch and process data from a **local instance** using InfraMock, there is nothing to do. The worker is
already configured to fetch data through it.

If you want to fetch and process data from **the live** `staging` or `production` clusters, you need to set
the `CLUSTER_ENV` environment variable to be either `staging` or `production`, respectively. The name of the SQS queue used
by the worker to subsribe to work is stored in the `WORKER_QUEUE` environment variable. You can set it to any active queue
that you know is currently being used, but it defaults to `development-queue.fifo`. This is created specifically and only
for local testing purposes.

In order to push work to the worker you can either:
* a) Submit work via the AWS console, using the SQS interface
* b) Submit work from the terminal using `aws-cli`

a) Just go to the AWS console, select the relevant queue and submit your work following the format described further in this section.

b) Make sure you have `aws-cli` installed. Then, put the work in the correct format, described further in this section in 
a file `payload.json`. Finally, use the `aws-cli` commands to send that payload to the desired SQS queue.
Here is an example:

    aws --endpoint-url=http://localhost:4566 sqs send-message --queue-url http://localhost:4566/queue/development-queue.fifo --message-body "$(< payload.json)" --message-group-id "$(date +%s)"

which will push the payload in `payload.json` to a local worker. You can do something similar if you are connected to the
live clusters by deleting the `--endpoint-url` argument and using the appropriate remote queue URL.

### Task formatting
To see what tasks are available for submitting what what is the format of each of them, go to the [API schema](https://github.com/biomage-ltd/api/blob/master/src/specs/api.yaml) and look at
`WorkRequest`. There, the schemas and parameters available and required are specified for all supported tasks.

### Sample tasks
Here are some examples:

* `GetEmbedding`:

    {
        "uuid": "509520fe-d329-437d-8752-b5868ad59425",
        "socketId": "Y1poEygzBfrDmIWpAAAA",
        "experimentId": "5e959f9c9f4b120771249001",
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
        "experimentId": "5e959f9c9f4b120771249001",
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
        "experimentId": "5e959f9c9f4b120771249001",
        "timeout": "2099-12-31 00:00:00",
        "body": {
            "name": "GeneExpression",
            "cellSets": ["louvain-0", "louvain-1", "louvain-2", "louvain-3", "louvain-4", "louvain-5", "louvain-6"],
            "genes": ["TGFB1", "CST3"]
        }
    }

* `PrepareExperiment`:

    {
        "uuid": "509520fe-d329-437d-8752-b5868ad59425",
        "socketId": "Y1poEygzBfrDmIWpAAAA",
        "experimentId": "1234",
        "timeout": "2099-12-31 00:00:00",
        "body": {
            "name": "PrepareExperiment",
            "sourceBucket": "biomage-source-originals",
            "sourceMatrixPath": "pbmc_count_matrices/hg19/"
        }
    }


### 3. Run the code
To run the worker locally, execute the following commands in a terminal:

    cd src/
    python work.py

Note that the worker will automatically switch itself off if it doesn't receive any tasks for 20 minutes. In this case, simply rerun it again and it will start, as normal.

## Run tests
Go to the src/ and run:

    CLUSTER_ENV="test" python -m pytest --cov=.