# worker-r

The Cellenics single cell data analysis tasks executor, written in R.

## Overview

The R part of the worker runs the data analysis tasks. It fullfills the following functions:

- Loads the processed `SeuratObject` for the experiment id assigned to the pod, reloading it if the associated file changes.
- Creates a local REST API specific to the assigned experiment id with [RestRserve](https://restrserve.org/index.html) that listens for requests from the Python worker.
- Runs the task function specified in the request from the Python worker.
- Returns the results back to the Python worker.

## Setup

The R worker (as of Jan 2022) uses R version **4.0.5**. To check your R version, you can run the following command in your terminal

    R --version

If you have a different version of R, please uninstall it, using the instructions in the [R manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Uninstalling-under-macOS). Once the package is uninstalled, download the required R version from the [CRAN hosted site](https://cran.r-project.org/bin/). For macOS, the required binaries are [here](https://cran.r-project.org/bin/macosx/base/). **Using the Homebrew distribution of R is not recommended.**

Some packages on MacOS X requires the Fortran compiler to compile. According to the [tooling page](https://cran.r-project.org/bin/macosx/tools/), R on Mac OSX requires GNU Fortran 8.2, which can be download from [gFortran for MacOS releases page](https://github.com/fxcoudert/gfortran-for-macOS/releases).

### Development with RStudio

The R worker is provided as a RStudio project, complete with a `renv`
definition. It might be useful to run things interactively in certain
development scenarios; to do so, you should have the correct R version installed (check `renv.lock` file for it). Open the `.Rproj` file with Rstudio and run in the R terminal:

    renv::restore()

In addition, dependencies for interactive development are specified in the DESCRIPTION file, under
the `Suggests` section. They, and their dependencies are ignored by renv (specified in the
[renv settings file](renv/settings.dcf) file).

They can be installed automatically by running the next block in an R console.

```{r}
    renv::install() # (without arguments)
```

### Development with Visual Studio Code

Aside from RStudio, development can also be done inside a development container that is automatically built, run, and managed by Visual Studio Code. You must have the [Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed, as specified in the main README.

To get started, open the r workspace:

    code r.code-workspace

You should be prompted to run the workspace inside a container. Accept this. Once you see the folder structure, the worker is running and you have access to the R worker's container.

To make sure everything works, try to access http://localhost:4000/health from your browser.This should give you a 200 HTTP response, which means the server is up.

### Running R worker Tests

_With R Studio_

Install the worker R package by running the following command in RStudio:

    install.packages('./')

Once the worker package is installed, run

    devtools::test()

from the RStudio terminal or press `cmd+shift+T`

_Via CLI_

- Go to worker/r folder
- Start an R session (enter `R` in the terminal)
- Run `devtools::test()`
