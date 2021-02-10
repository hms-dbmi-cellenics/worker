#!/usr/bin/env Rscript

# Install R packages or fail with error.
#
# Arguments:
#   - First argument (optional) can be one of:
#       1. repo URL
#       2. "github" if installing from GitHub repo (requires that package 'devtools' is
#          already installed)
#       3. full URL of package from which to install from source; if used, provide package
#          name in second argument (e.g. 'curl')
#     If this argument is omitted, the default repo https://cran.rstudio.com/ is used.
#   - Remaining arguments are either:
#       1. one or more R package names, or
#       2. if installing from GitHub, the path containing username and repo name, e.g.
#          'timelyportfolio/htmlwidgets_spin', optionally followed by the package name (if
#          it differs from the GitHub repo name, e.g. 'spin').

arg_list <- commandArgs(trailingOnly = TRUE)

if (length(arg_list) < 1) {
    print("ERROR: Too few arguments.")
    quit(status = 1, save = "no")
}

if (arg_list[1] == "github" || grepl("^https?://", arg_list[1], perl = TRUE)) {
    if (length(arg_list) == 1) {
        print("ERROR: No package name provided.")
        quit(status = 1, save = "no")
    }
    repo <- arg_list[1]
    packages <- arg_list[-1]
} else {
    repo <- "https://cran.rstudio.com/"
    packages <- arg_list
}

for (i in seq_along(packages)) {
    p <- packages[i]

    start_time <- Sys.time()
    if (grepl("^https?://[A-Za-z0-9.-]+/.+\\.tar\\.gz$", repo, perl = TRUE)) {
        # If 'repo' is URL with path after domain name, treat it as full path to a package
        # to be installed from source.
        install.packages(repo, repo = NULL, type = "source")
    } else if (repo == "github") {
        # Install from GitHub.
        github_path <- p
        elems <- strsplit(github_path, "/")
        if (lengths(elems) != 2) {
            print("ERROR: Invalid GitHub path.")
            quit(status = 1, save = "no")
        }
        username <- elems[[1]][1]
        github_repo_name <- elems[[1]][2]
        if (!is.na(packages[i + 1])) {
            # Optional additional argument was given specifying the R package name.
            p <- packages[i + 1]
        } else {
            # Assume R package name is the same as GitHub repo name.
            p <- github_repo_name
        }

        library(devtools)
        install_github(github_path)
    } else {
        # Install from R package repository.
        install.packages(p, dependencies = TRUE, repos = repo)
    }
    end_time <- Sys.time()

    if (!library(p, character.only = TRUE, logical.return = TRUE)) {
        quit(status = 1, save = "no")
    } else {
        cat(paste0("Time to install ", p, ":\n"))
        print(end_time - start_time)
    }

    if (repo == "github") {
        break
    }
}