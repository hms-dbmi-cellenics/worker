worker-python
=============

The Cellenics single cell analysis tasks wrapper, written in Python.

Overview
--------
The Python part of the worker is the main entry point for data analysis tasks. It fullfills the following functions:
- Receives tasks from the API by listening to an SQS queue associated with the relevant experiment ID.
- Prepares the received task for the R part of the worker. This is task dependent and may include additional task validation, cleaning and formatting the data or fetching additional data from AWS services.
- Forwards the task to the R part of the worker, using the local network.
- Receives the results back from the R part of the worker, uploads the data to S3 and sends a notification to Redis that the task has been computed, using [socket io API](https://pypi.org/project/socket.io-emitter/).

Running the worker
------------------

### Start the worker

See the main README for instructions on how to run the workers in Docker.

### Process tasks
Tasks are automatically processed when they are received from the SQS queue specified.

For **local development**, make sure you have [InfraMock](https://github.com/hms-dbmi-cellenics/inframock)
running alongside the [ui](https://github.com/hms-dbmi-cellenics/ui) and [api](https://github.com/hms-dbmi-cellenics/api)
projects. Refer to their respective documentations on how to run them locally. Once all of these are running,
tasks should automatically be submitted and processed when you perform actions on the `ui`. There is nothing else to do.


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
Typically, you will be able to run `make test` in the outer directory to execute the unit tests. If you need
to run the tests manually, for instance if you only want to run a specific test, do the following:

Either:
- open a terminal in VS Code while the workspace is open in the container (see the previous section) 
- shell into the docker container manually (`docker exec -it biomage-worker-python bash`)

Then run these commands:

    export CLUSTER_ENV="development"
    python -m pytest --cov=. --cov-report=term-missing

You can append specific test files (e.g. `tests/tasks/test_cluster_cells.py`) to the pytest command to only
run those tests.

### Task formatting
Task definitions are stored in the `api` project as an OpenAPI schema.
You can find this [here](https://github.com/hms-dbmi-cellenics/api/blob/master/src/specs/api.yaml).

Download the schema and open it using Stoplight Studio. Looking into `WorkRequest` should give you the schemas and parameters for all supported tasks.