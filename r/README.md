worker-r
========

The Cellenics single cell analysis tasks executor, written in R.

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
