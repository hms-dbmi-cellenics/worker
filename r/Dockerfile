# Create builder step
FROM rocker/r-ver:4.4.0 AS builder
WORKDIR /src/worker

# install required debian packages to install R packages
COPY setup/install_debian_packages.sh .
COPY setup/sysdeps_build_debian.txt .
RUN cat sysdeps_build_debian.txt | xargs ./install_debian_packages.sh

# need GITHUB_PAT to authenticate github installations
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -q -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# add renv library to .libPaths
# this is path with R package folders when renv deactivated or not used
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install renv to install required R packages
RUN R -q -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -q -e "remotes::install_github('rstudio/renv@v1.0.7')" && \
    R -q -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# fast pre-restore with pkgbuild
COPY ./renv.lock.init ./renv.lock
COPY setup/restore_fast.R .
RUN Rscript restore_fast.R

# restore anything pkgbuild missed
COPY setup/restore_renv.R .
RUN Rscript restore_renv.R

# use renv::snapshot() while R dependency updates are quick to build
COPY ./renv.lock .
RUN Rscript restore_renv.R

# deactivate so that above .libPaths used
RUN R -q -e "renv::deactivate()"
RUN R -q -e "remotes::install_github('bnprks/BPCells/r@c293a0a')"

# determine system run-time deps
COPY setup/get_sysdeps_run.R .
RUN Rscript get_sysdeps_run.R

# remove identified packages with potentially conflicting licenses
COPY setup/create_ghost_packages.R .
RUN Rscript create_ghost_packages.R

# check for any potentially problematic licenses
COPY setup/check_package_licenses.R .
RUN Rscript check_package_licenses.R

# ---------------------------------------------------
# COMMON MINIMAL BUILD
# ---------------------------------------------------
FROM rocker/r-ver:4.4.0 AS common
WORKDIR /src/worker

# get source code and R packages
COPY --from=builder /src /src

# add path where renv installed R packages to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install runtime system deps
# python3-dev prevents 'reticulate can only bind to copies of Python built with --enable-shared'
# cleanup setup files
RUN echo "python3-pip python3-venv python3-dev zstd" >> sysdeps_run.txt && \
    cat sysdeps_run.txt | xargs ./install_debian_packages.sh && \
    rm -rf *

ENV WORKON_HOME=/src/.virtualenvs
RUN R -q -e "reticulate::virtualenv_create('r-reticulate', python='$(which python3)')" && \
    R -q -e "reticulate::virtualenv_install('r-reticulate', c('setuptools==68.2.2', 'umap-learn==0.5.3', 'numpy<2'), pip_options='--no-cache-dir --upgrade')"

# ---------------------------------------------------
# CASSIA (LLM-based cell type annotation)
# Installed at the very end so it doesn't invalidate the expensive renv /
# BPCells layers above. CASSIA is a reticulate wrapper around a bundled Python
# backend, so it needs its python deps AND the R package (from GitHub below).
# The python deps go into the same r-reticulate venv the worker already binds
# (reticulate binds a single Python per process); numpy is kept <2 for umap-learn.
# ---------------------------------------------------
# boto3 + aws-bedrock-token-generator let the worker mint a short-term Bedrock
# bearer token from its IAM role (IRSA) at runtime, so no static API key is
# needed. openai is the client CASSIA's custom provider uses to call Bedrock's
# OpenAI-compatible endpoint.
RUN R -q -e "reticulate::virtualenv_install('r-reticulate', c('numpy<2', 'pandas', 'openai', 'anthropic', 'requests', 'matplotlib', 'seaborn', 'mygene', 'boto3', 'aws-bedrock-token-generator'), pip_options='--no-cache-dir')"

# CASSIA R package (reticulate wrapper + bundled python backend), pinned commit.
# GITHUB_PAT avoids the unauthenticated GitHub API rate limit during install.
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT
RUN R -q -e "remotes::install_github('ElliotXie/CASSIA/CASSIA_R@7a77067f2af4ac89e68e8424a80ccd5e3dc8c334')"

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
RUN pip install -U jedi PyYAML watchdog[watchmedo]

# add R package and runner
ADD R ./R
ADD tests ./tests
COPY DESCRIPTION NAMESPACE work.R renv.lock ./

CMD ["Rscript", "work.R"]
