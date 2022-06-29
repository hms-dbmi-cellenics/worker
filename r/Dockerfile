# Create builder step
FROM rocker/r-ver:4.2.0 AS builder
WORKDIR /src/worker

# install required debian packages to install R packages
COPY setup/install_debian_packages.sh .
COPY setup/sysdeps_build_debian.txt .
RUN cat sysdeps_build_debian.txt | xargs ./install_debian_packages.sh

# need GITHUB_PAT to authenticate github installations
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install renv to install required R packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.15.5')" && \
    R -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# an initial lockfile is used to avoid frequent re-installs
# use renv::snapshot(lockfile='renv.lock.init') if R dependency updates become slow to build
# delete renv cache
# strip debug from shared libraries
# see http://dirk.eddelbuettel.com/blog/2017/08/20/#010_stripping_shared_libraries
COPY ./renv.lock.init .
RUN R -e "renv::restore(lockfile='renv.lock.init', library = '$RENV_LIB')" && \
    R -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so

RUN R -e "renv::deactivate()"

# install miniconda and python umap-learn for RunUMAP
# clean conda
ENV RETICULATE_MINICONDA_PATH=/src/r-miniconda
RUN R -e "reticulate::install_miniconda()" && \
    R -e "reticulate::conda_install(packages = 'umap-learn=0.5.3', python_version='3.10')" && \
    CONDA_PATH=$(echo "cat(reticulate::conda_binary())" | Rscript -) && \
    $CONDA_PATH clean --force-pkgs-dirs -y

# use renv::snapshot() while R dependency updates are quick to build
COPY ./renv.lock .
RUN R -e "renv::restore(lockfile='renv.lock', library = '$RENV_LIB', clean = TRUE)" && \
    R -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so

# determine system run-time deps
COPY setup/get_sysdeps_run.R .
RUN Rscript get_sysdeps_run.R

# ---------------------------------------------------
# COMMON MINIMAL BUILD
# ---------------------------------------------------
FROM rocker/r-ver:4.2.0 AS common
WORKDIR /src/worker
ENV RETICULATE_MINICONDA_PATH=/src/r-miniconda

# get source code and R packages
COPY --from=builder /src /src

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install runtime system deps
# cleanup setup files
RUN echo "python3-pip" >> sysdeps_run.txt && \
    cat sysdeps_run.txt | xargs ./install_debian_packages.sh && \
    rm -rf *

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM common AS prod

# add R package and runner
ADD R ./R
ADD tests ./tests
COPY DESCRIPTION NAMESPACE work.R ./

# start app
ENTRYPOINT ["bash", "/var/lib/watchfile/entrypoint.sh"]
CMD ["Rscript", "work.R"]

# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM common AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

# add R package and runner
ADD R ./R
ADD tests ./tests
COPY DESCRIPTION NAMESPACE work.R ./

CMD watchmedo auto-restart --directory=. --pattern='*.R' --recursive -- Rscript work.R
