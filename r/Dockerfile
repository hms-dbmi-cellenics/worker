# Create builder step
FROM r-base:4.0.5 AS builder
WORKDIR /src

# install required debian packages to make R packages work
COPY setup ./setup
COPY requirements_debian.txt .
RUN cat requirements_debian.txt | xargs ./setup/install_debian_packages.sh
RUN apt-get update && apt-get install -y procps

# need GITHUB_PAT to authenticate github installations
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# install required R packages for the worker.
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.13.2')"

# an initial lockfile is used to avoid frequent re-installs
# use renv::snapshot(lockfile='renv.lock.init') if R dependency updates become slow to build
COPY renv.lock.init .
RUN R -e "renv::init(bare = TRUE); options(renv.consent = TRUE); renv::restore(lockfile='renv.lock.init')"

# use renv::snapshot() while R dependency updates are quick to build
COPY renv.lock .
RUN R -e "options(renv.consent = TRUE); renv::restore(lockfile='renv.lock')"

RUN R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor','Matrix.utils'))"

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('cole-trapnell-lab/leidenbase')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3')"

# install miniconda and python umap-learn for RunUMAP
ENV RETICULATE_MINICONDA_PATH "/src/renv/python/r-reticulate"
RUN R -e "reticulate::install_miniconda()"
RUN R -e "reticulate::conda_install(Sys.getenv('RETICULATE_MINICONDA_PATH'), packages='umap-learn')"
RUN R -e "install.packages('covr')"
RUN R -e "install.packages('xml2')"


RUN apt-get install -y libudunits2-dev
RUN apt-get install -y libgdal-dev

RUN R -e "install.packages('sf')"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('cole-trapnell-lab/leidenbase')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3')"
RUN R -e "install.packages('R.utils')"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"


# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# copy R source files and runner
ADD R ./R
COPY work.R .

# start app
ENTRYPOINT ["bash", "/var/lib/watchfile/entrypoint.sh"]
CMD ["Rscript", "work.R"]

# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

# copy R source files and runner
ADD R ./R
COPY work.R .

CMD watchmedo auto-restart --directory=. --pattern='*.R' --recursive -- Rscript work.R