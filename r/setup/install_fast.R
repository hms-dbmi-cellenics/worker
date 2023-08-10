renv::install('pkgdepends')
lockfile <- renv:::renv_lockfile_read('renv.lock')

# get records to install
records <- renv:::renv_lockfile_records(lockfile)

# convert into specs compatible with pak, and install
remotes <- renv:::map_chr(records, renv:::renv_record_format_remote)

urls <- lapply(records, `[[`, 'RemoteUrl')
has.url <- !sapply(urls, is.null)
remotes[has.url] <- paste0('url::', unlist(urls))

# to get dependencies in order that they need to be installs
get_ordered_recursive_dependencies <- function(lockfile, target_package) {
  ordered_dependencies <- character(0)
  visited_packages <- character(0)

  find_dependencies <- function(pkg_name) {
    if (pkg_name %in% visited_packages || !(pkg_name %in% names(lockfile$Packages))) {
      return()
    }

    visited_packages <<- c(visited_packages, pkg_name)

    pkg_info <- lockfile$Packages[[pkg_name]]
    pkg_requirements <- pkg_info$Requirements

    for (req_pkg in pkg_requirements) {
      find_dependencies(req_pkg)
    }

    if (!(pkg_name %in% ordered_dependencies)) {
      ordered_dependencies <<- c(ordered_dependencies, pkg_name)
    }
  }

  find_dependencies(target_package)

  return(ordered_dependencies)
}


# try all at once
# invalid version specification 'NA' was from Matrix.utils
# too resolve: see what isn't available
# avail <- available.packages()
# not.avail <- names(records)[!names(records) %in% row.names(avail)]
# renv::install(url/to/archive)
# renv::snapshot()
prop <- pkgdepends::new_pkg_installation_proposal(
  remotes, config = list(library = .libPaths()[[1]], dependencies = FALSE), policy = 'lazy')

prop$download()
plan <- prop$get_install_plan()

fix_plan <- function(plan, lockfile) {
  deps <- list()
  for (i in 1:nrow(plan)) {
    package <- plan$package[i]
    depsi <- lockfile$Packages[[package]]$Requirements
    deps[[i]] <- depsi
  }
  plan$dependencies <- deps
  return(plan)
}
plan <- fix_plan(plan, lockfile)

pkgdepends:::install_package_plan(plan, num_workers = pkgdepends:::get_num_workers())

install_package_plan <- function (plan, lib = .libPaths()[[1]], num_workers = 1, cache = NULL)
{
  start <- Sys.time()
  cli::ansi_hide_cursor()
  on.exit(cli::ansi_show_cursor())
  cli::cli_div(theme = list(.timestamp = list(color = "darkgrey",
                                              before = "(", after = ")")))
  required_columns <- c("type", "binary", "dependencies",
                        "file", "needscompilation", "package")
  pkgdepends:::assert_that(inherits(plan, "data.frame"), all(required_columns %in%
                                                  colnames(plan)), is_string(lib), is_count(num_workers,
                                                                                            min = 1L))
  if (!"vignettes" %in% colnames(plan))
    plan$vignettes <- FALSE
  if (!"metadata" %in% colnames(plan)) {
    plan$metadata <- replicate(nrow(plan), character(),
                               simplify = FALSE)
  }
  if (!"packaged" %in% colnames(plan))
    plan$packaged <- TRUE
  plan <- add_recursive_dependencies(plan)
  config <- list(lib = lib, num_workers = num_workers, show_time = tolower(Sys.getenv("PKG_OMIT_TIMES")) !=
                   "true")
  state <- make_start_state(plan, config)
  state$cache <- cache
  state$progress <- create_progress_bar(state)
  on.exit(done_progress_bar(state), add = TRUE)
  withCallingHandlers({
    for (i in seq_len(state$config$num_workers)) {
      task <- select_next_task(state)
      state <- start_task(state, task)
    }
    repeat {
      if (are_we_done(state))
        break
      update_progress_bar(state)
      events <- poll_workers(state)
      state <- handle_events(state, events)
      task <- select_next_task(state)
      state <- start_task(state, task)
    }
  }, error = function(e) kill_all_processes(state))
  create_install_result(state)
}



# try grouped
for (i in seq_along(remotes)) {
  pkg <- remotes[i]
  target_package <- names(pkg)
  deps <- get_ordered_recursive_dependencies(lockfile, target_package)
  pkgs <- remotes[deps]

  cat('DOWNLOADING:', pkgs, '...\n')
  # perform installation
  prop <- pkgdepends::new_pkg_installation_proposal(
    pkgs, config = list(library = .libPaths()[[1]], dependencies = FALSE), policy = 'lazy')

  prop$download()
}

ordered_deps <- get_ordered_deps(lockfile)
ordered_deps <- remotes[ordered_deps]

for(pkg in ordered_deps) {
  cat('DOWNLOADING:', pkg, '...\n')
  # perform installation
  prop <- pkgdepends::new_pkg_installation_proposal(
    pkg, config = list(library = .libPaths()[[1]], dependencies = FALSE))

  try(prop$download())
}
#


blah <- renv:::renv_retrieve_impl(pkg)

type <- 'source'
repos <- lockfile$R$Repositories

blah <- renv:::renv_available_packages_query(type, repos)
