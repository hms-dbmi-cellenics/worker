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
    R -e "remotes::install_github('rstudio/renv@0.14.0')"

# an initial lockfile is used to avoid frequent re-installs
# use renv::snapshot(lockfile='renv.lock.init') if R dependency updates become slow to build
COPY renv.lock.init .
RUN R -e "renv::activate(); renv::restore(lockfile='renv.lock.init')"

# use renv::snapshot() while R dependency updates are quick to build
COPY renv.lock .
RUN R -e "renv::restore(lockfile='renv.lock')"

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

CMD watchmedo auto-restart --directory=./R --pattern=* --recursive -- Rscript work.R
