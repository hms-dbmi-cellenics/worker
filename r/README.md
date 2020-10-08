worker-r
========

The single-cell pipeline work executor for R.

Setup
-----

Open the r workspace:

    code r.code-workspace

Development is done inside a development container that is automatically built,
run, and managed by Visual Studio Code. You do not need R, R Studio, or a Python
virtual environment to be installed locally.

You must have the following the [Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed, as specified
in the main README.

You should be prompted to run the workspace inside a container. Accept this. Once
you see the folder structure, the worker is running and you have access to the
R worker's container.

To make sure everything works, try to access http://localhost:4000 from your browser.
The development container automatically forwards this port while you are attached
to the container.