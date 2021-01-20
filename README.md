# worker

The single-cell pipeline work executor.

## Running locally

While in the `worker/` root folder on the host, you can use `docker-compose` as you would normally.

For example, to run the r and python containers, you can do:

    docker-compose up --build

On Linux, this instead needs to be:

    docker-compose -f docker-compose.linux-dev.yaml up --build

Note that during the first time, the build can take up to 40-50 minutes to complete.
If you get an error, see the `Troubleshoooting` section for help.

To get a development log stream of both containers running, you can use:

    docker-compose logs -f

To shut down the development containers, you can use:

    docker-compose kill

## Development

### Prerequisites

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
container, try running `docker-compose up --build` to see where exactly the build breaks.
Please check `Troubleshooting` section that lists commonly occuring problems.

The root directories of each of the workspaces are dynamically linked to `/r` and `/python`
respectively. The terminals spawn terminals within the containers, as expected.

These development environments should be pre-configured with the same requirements as the
produciton instances, as well as the necessary VS Code extensions required to debug and
lint code.

## More details

For more details on the individual runners, check out the README files in their respective directories.

## Troubleshooting

1.  Errors saying `... unsupported option: 'target'` after running `docker-compose up --build`.

    This is most likely a problem with the docker-compose version. Simply re-install it:

         pip3 uninstall docker-compose
         pip3 install -U docker-compose

2.  `docker-compose up --build` fails due to rate limit errors.

    To fix this one, make sure you create a personal access token in your Github account and
    add it as an environment variable, called `GITHUB_PAT`:

        1. Go to https://github.com/settings/tokens, create a new token. The token should be read only.
        2. Set `GITHUB_PAT` to equal to the value of the token in a terminal.

3.  `Error: Failed to install 'unknown package' from Github: Timeout was reached: [api.github.com] Resolving timed out after 10000 milliseconds`

    This error is due to a bug in DNS resolution of Alpine-based containers running on early releases of Docker Desktop for Mac version 3.

    To fix this, you can download and use a previous version of Docker (e.g. 2.5.0.1) from https://docs.docker.com/docker-for-mac/release-notes/
