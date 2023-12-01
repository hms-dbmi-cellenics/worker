[![codecov](https://codecov.io/gh/hms-dbmi-cellenics/worker/branch/master/graph/badge.svg?token=3PHqr61GpH)](https://codecov.io/gh/hms-dbmi-cellenics/worker)
worker
======

The Cellenics data analysis tasks executor. It consists of two containers: a Python container and an R container. For more details on the individual containers, check out the README files in their respective directories.

The Python part of the worker is a wrapper around the R part: it receives tasks from the API, parses them, sends them to the R part for computation, then formats the results, uploads them to S3 and sends a notification to the API via Redis-backed socket connection that they are ready.

The R part of the worker computes single cell analysis tasks on a pre-processed Seurat rds object, loaded into memory from S3. The R part of the worker can communicate only with the Python part of the worker.

More specific details about the Python or the R part of the worker can be found in the README file in the respective folder (python/ or r/).

## Deployment

The worker is deployed as a Helm chart to an AWS-managed Kubernetes cluster and runs on a Fargate-managed node. The Helm chart template for the worker is located in `chart-infra/` folder.

 The deployment of the worker is handled by the cluster Helm operator and the [worker Github Actions workflow](https://github.com/hms-dbmi-cellenics/worker/blob/master/.github/workflows/ci.yaml).

During a deployment, the worker Github Actions workflow does the following:
- It pushes new worker images to ECR.
- Adds deployment-specific configurations to the worker Helm chart. Pushes those deployment-specific configuration changes in [releases/](https://github.com/hms-dbmi-cellenics/iac/tree/master/releases) folder in iac, under the relevant environment.

## Development

### Prerequisites

#### Docker resource allocation

Make sure that sufficient resources are allocated in Docker to be able to
compile everything. 10gb RAM and 20gb disk image size should be more than
enough. Be mindful when allocating RAM, too much and you could end up freezing
your computer.

#### Github Access Token

Create a [Github Access
Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
and export it as an environment variable before building. You should consider
adding it to your `.bashrc` (or `.zshrc`) for convenience.

``` shell
export GITHUB_API_TOKEN=<your-token>
```

### Running tests

#### python worker
Assuming the containers are running, you can execute the (pytest) unit tests using:

    make test

See [here](python/README.md#tests) for more information about the tests.

To shut down the development containers, you can use:

    make kill

#### R worker

Refer to the [R-worker README](r/README.md#running-r-worker-tests).

### Remote - Containers
Development is done inside a development container that is automatically built,
run, and managed by Visual Studio Code. You do not need R, R Studio, or a Python
virtual environment to be installed locally.

As such, you must have the
[Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
extension installed. Make sure you restart VS Code after installing to make sure it
loads successfully. You should see a green icon in the leftmost part of the status bar,
which indicates that the remote container plugin has been installed.

### Setup

To open the R workspace, you can type `code r/r.workspace` while in the terminal inside
VS Code.

Similarly, to open the Python workspace, you can type `code python/python.code-workspace`.

You should be prompted to run the workspace inside a container. Accept this. Once
you see the folder structure, the worker is running and you have access to the
R worker's container. If you get an error after trying to run the workspace inside a
container, try running `make build` to see where exactly the build breaks.
Please check `Troubleshooting` section that lists commonly occuring problems.

The root directories of each of the workspaces are dynamically linked to `/r` and `/python`
respectively. The terminals spawn terminals within the containers, as expected.

These development environments should be pre-configured with the same requirements as the
produciton instances, as well as the necessary VS Code extensions required to debug and
lint code.

## Running locally
To run the worker locally you need to:
1. Build the docker image locally (1.1) or Download the docker image directly from ECR (1.2)
2. Run it, passing the id of the processed experiment that you want to use the worker with.

### 1. Obtaining worker image
#### 1.1. Building an image locally
While in the `worker/` root folder on the host, you can use `make build`.

To build and run the R and python containers, you can do:

    make build

Note that during the first time, the build can take up to 40-50 minutes to complete.
If you get an error, see the `Troubleshoooting` section for help.

To get a development log stream of both containers running, you can use:

    make logs

#### 1.2 Downloading the image from ECR
While in the `worker/` root folder on the host, you can use:
```
AWS_REGION={Your AWS region} AWS_ACCOUNT_ID={Your account ID} LATEST_TAG={Latest release tag} make download-image make download-image.
```

### Step 2. Running the worker
Before running the worker, you need to have a folder named with the experiment id that you want to load. The folder should be saved under `worker/data` and it has to contain:
 - Processed rds object file, called `r.rds`
 - A json file of the cell_sets for that experiment, called `cell_sets.json`

Here is an example folder structure for experiment id `1234`:

```bash
data
├── 1234
│   ├── cell_sets.json
│   └── r.rds
```

You can obtain this folder structure if you do either of the following:
- Start the rest of the Cellenics components locally, then upload samples from Data Management and launch analysis. A prerequisite for this step is configuring the pipeline to run locally.

OR

- Download the r.rds object and the cell_sets.json file for an already processed experiment from S3 and then manually add them in `worker/data` under a new folder named with the experiment id.

You can have one or more experiments under `worker/data`.

To run the worker with the experiment id of your choice, do the following:

**If the image is build locally**
In a terminal, while in the `worker/` root folder, type the following:

    EXPERIMENT_ID=1234 make run

where `1234` is the experiment id of your choice.

**If the image is downloaded from ECR**
In a terminal, while in the `worker/` root folder, type the following:

    EXPERIMENT_ID=1234 make run-downloaded

where `1234` is the experiment id of your choice.


## Troubleshooting

1.  Errors saying `... unsupported option: 'target'` after running `make build`.

    This is most likely a problem with the docker-compose version (which is used by `make` to build the images). Simply re-install it:

         pip3 uninstall docker-compose
         pip3 install -U docker-compose

2.  `make build` fails due to rate limit errors.

    To fix this one, make sure you create a personal access token in your Github account and
    add it as an environment variable, called `GITHUB_PAT`:

        1. Go to https://github.com/settings/tokens, create a new token. The token should be read only.
        2. Set `GITHUB_PAT` to equal to the value of the token in a terminal.

3.  `Error: Failed to install 'unknown package' from Github: Timeout was reached: [api.github.com] Resolving timed out after 10000 milliseconds`

    This error is due to a bug in DNS resolution of Alpine-based containers running on early releases of Docker Desktop for Mac version 3.

    To fix this, you can download and use a previous version of Docker (e.g. 2.5.0.1) from https://docs.docker.com/docker-for-mac/release-notes/

4.  When working locally only: `Keyboard interrupt` shows up non-stop, impeding the worker to ever start up.

    We believe this error is related to the hot reload mechanism and probably having some temporary files that, on changing, trigger the hot reload non-stop. The workaround to be able to continue working is to disable hot-reload (so, for changes to go through, you'll need to restart the worker manually). Just performing [these changes](https://github.com/hms-dbmi-cellenics/worker/pull/76/files) in your local Dockerfiles should be enough

5.  Error when attempting to start the worker saying something like:
`botocore.exceptions.EndpointConnectionError: Could not connect to the endpoint URL: "http://host.docker.internal:4566/biomage-source-development?...`

First, check inframock is running. If it isn't, start it and try again. Otherwise, see below.

##### For Linux users

**Note** this error should already been handled by the Makefile builds. If you encounter it while
using `make run`, report it in the Slack channel #engineering.

Go to `docker-compose.yaml`
In the `python` and `r` entries add at the end:

```
extra_hosts:
      - "host.docker.internal:host-gateway"
```

IMPORTANT: Don't include this in a PR, because it will break stuff on macOS.

### Debugging locally

**TLDR:** Save anything in /debug in the container and it will be available at `$(pwd)/data/debug`.

To save the `req` argument to a worker function, specify DEBUG_STEP. DEBUG_STEP can be either `all` (will save `req` from any task) or the basename of a [path in work.R](r/work.R#L42) and will hot-reload if changed in work.R. It can also be set on initial run:

```bash
# e.g. DEBUG_STEP=getClusters
DEBUG_STEP=task_name make build && make run
```

When a worker function is run, it will save the `req` and `data` objects used by the specified `task_name` in ./data/debug. You will see a prompt to
read these into your R environment:

```R
# clicking these files in RStudio does this for you
req <- readRDS('./data/debug/{experiment_id}_{task_name}_req.rds')
data <- readRDS('./data/debug/{experiment_id}_data.rds')
```
