# worker

The single-cell pipeline work executor.

## Running locally

While in the `worker/` root folder on the host, you can use `make build`.

For example to build and run the r and python containers, you can do:

    make build && make run


Note that during the first time, the build can take up to 40-50 minutes to complete.
If you get an error, see the `Troubleshoooting` section for help.

To get a development log stream of both containers running, you can use:

    make logs

To shut down the development containers, you can use:

    make kill

## Development

### Prerequisites

#### Remote - Containers
Development is done inside a development container that is automatically built,
run, and managed by Visual Studio Code. You do not need R, R Studio, or a Python
virtual environment to be installed locally.

As such, you must have the
[Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
extension installed. Make sure you restart VS Code after installing to make sure it
loads successfully. You should see a green icon in the leftmost part of the status bar,
which indicates that the remote container plugin has been installed.

#### Git LFS
File(s) under `data/test` are downloaded by [inframock](https://github.com/biomage-ltd/inframock), uploaded to mock S3 and used by the workers. As some of these files are over Github's file size limit (100 MB), they are stored using [Git LFS](https://git-lfs.github.com/). Follow the installation instructions on their website to setup Git LFS locally.

Once you have installed Git LFS, you can open the worker root directory in a terminal and run 

```
git lfs install
```

If Git LFS is installed successfully, it should print

```
Updated git hooks.
Git LFS initialized.
```

You can see the list of files tracked by Git LFS in `.gitattributes`.

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

## More details

For more details on the individual runners, check out the README files in their respective directories.

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
    
### For linux users

1. Error when attempting to start the worker saying something like:
`botocore.exceptions.EndpointConnectionError: Could not connect to the endpoint URL: "http://host.docker.internal:4566/biomage-source-development?...`

**Note** this error should already been handled by the Makefilebuilds . If you encounter it while using `make build` report in on slack channel #engineering.

Go to `docker-compose.yaml`
In the python and r entries add at the end:

```
extra_hosts:
      - "host.docker.internal:host-gateway"
```

IMPORTANT: Don't include this in a PR, because it will break stuff on macOS.



## Debugging locally

To save the `req` argument to a worker function, specify DEBUG_STEP. DEBUG_STEP can be either `all` (will save `req` from any task) or the basename of a [path in work.r](r/src/work.r#L88) and will hot-reload if changed at the top of work.r. It can also be set on initial run:

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
