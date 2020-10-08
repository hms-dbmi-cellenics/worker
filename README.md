worker
======

The single-cell pipeline work executor.

Setup
-----

Development is done inside a development container that is automatically built,
run, and managed by Visual Studio Code. You do not need R, R Studio, or a Python
virtual environment to be installed locally.

As such, you must have the
[Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
extension installed. Make sure you restart VS Code after installing to make sure it
loads successfully. You should see a green icon in the leftmost part of the status bar,
which indicates that the remote container plugin has been installed.

Development
-----------

To open the R workspace, you can type `code r/r.workspace` while in the terminal inside
VS Code.

Similarly, to open the Python workspace, you can type `code python/python.code-workspace`.

You should be prompted to run the workspace inside a container. Accept this. Once
you see the folder structure, the worker is running and you have access to the
R worker's container. In the bottom right corner you should see the name of the container
VS Code is running in.

The root directories of each of the workspaces are dynamically linked to `/r` and `/python`
respectively. The terminals spawn terminals within the containers, as expected.

These development environments should be pre-configured with the same requirements as the
produciton instances, as well as the necessary VS Code extensions required to debug and
lint code.

Managing the containers
-----------------------

While in the `worker/` root folder on the host, you can use `docker-compose` as you would normally.

For example, to get a development log stream of both containers running, you can use:

    docker-compose logs -f

To shut down the development containers, you can use:

    docker-compose kill

More details
------------

For more details on the individual runners, check out the README files in their respective directories.