This repository contains 

# Overview
The purpose of the worker is to carry out data analysis tasks, for example computing an embedding. Data analysis tasks get received from a SQS queue and results will get submitted to a SNS topic the API is subscribed to.


# Build and Deploy step
The deployment step for this repository consits of deploying all resources needed for work to be carried out. Note that the deployment step in this repository does not start the work. Work gets started dynamically upon a user request by the API: https://gitlab.com/biomage/api.

CI deployment needs four variables that need to be created on the CI/CD side:

* `GITLAB_DEPLOY_USER` and `GITLAB_DEPLOY_PASSWORD`. These need to specify a valid GitLab deploy key that has read rights to the Docker registry. This secret is used by dynamically generated jobs to automatically be able to pull the worker image.
* `K8S_SECRET_AWS_ACCESS_KEY_ID`, `K8S_SERET_AWS_SECRET_ACCESS_KEY`, `K8S_SECRET_AWS_DEFAULT_REGION`, which specify the IAM account and default region to be used with this worker.

# Development

## Run locally

### 1. Set up environment
To run this code locally, first make sure you start a python virtual environment, activate it and install all the requirements:

        python3.7 -m venv
        source venv/bin/activate
        pip install -r requirements.txt

The code in the worker requires an access to AWS resources, so make sure you have aws cli installed and your machine has access to aws.

### 2. Submit a task
The next step is to send a task to the worker. To do that, you have to submit the desired task to the SQS queue it is subscribed to. By default, the name of the SQS queue used by the worker is stored in the environment variable `WORKER_QUEUE`. If `WORKER_QUEUE` is not defined, the default queue that will be polled is called `test-queue`.

To submit a task to `test-queue`

The "test-queue" is created specifically and only for local testing purposes. Go to the AWS console, under queues and select it (if it doesn't exist, create one using the user interface). 

To submit a GetEmbedding task, you can paste this in the SQS:

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

To submit a ListGenes task, you can paste this in the SQS:
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

it will be picked by the worker, the task will be computed and the results sent back via SNS.

Where `count_matrix` is the S3 bucket and key of the anndata file you want processed.

### 3. Run the code
After you have submitted a task, run:

    python src/worker.py
