worker-r
========

The Cellenics single cell data analysis tasks executor, written in R.

Overview
--------
The R part of the worker runs the data analysis tasks. It fullfills the following functions:
- Loads the processed `SeuratObject` for the experiment id assigned to the pod, reloading it if the associated file changes.
- Creates a local REST API specific to the assigned experiment id with [RestRserve](https://restrserve.org/index.html) that listens for requests from the Python worker.
- Runs the task function specified in the request from the Python worker.
- Returns the results back to the Python worker.
Setup
-----

Open the r workspace:

    code r.code-workspace

Development is done inside a development container that is automatically built,
run, and managed by Visual Studio Code. You do not need R, R Studio, or a Python
virtual environment to be installed locally.

You must have the[Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed, as specified
in the main README.

You should be prompted to run the workspace inside a container. Accept this. Once
you see the folder structure, the worker is running and you have access to the
R worker's container.

To make sure everything works, try to access http://localhost:4000/health from your browser.
This should give you a 200 HTTP response, which means the server is up.


### Running R worker interactively in Rstudio

The R worker is provided as a Rstudio project, complete with a `renv`
definition. It might be useful to run things interactively in certain
development scenarios; to do so, you should have the correct R version installed
(check `renv.lock.init` file for it). Open the `.Rproj` file with Rstudio and
run in the R terminal:

``` R
renv::restore()
```

### Running R worker Tests

In Rstudio:
devtools::test() or cmd+shift+T

Outside Rstudio:

Go to worker/r folder
Start an R session (R in the terminal)
Run devtools::test()
