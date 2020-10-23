# Create builder step
# pull official base image and use it as builder
FROM r-base:4.0.2 AS builder

WORKDIR /src

# install required debian packages to make R packages work
COPY setup ./setup
COPY requirements_debian.txt .
RUN cat requirements_debian.txt | xargs ./setup/install_debian_packages.sh

# install required R packages for the worker.
# these are sectioned as seen here so we can cache compiled code
# NOTE: add new packages to the END of this list. this makes sure
# previously installed packages don't get recompiled for no reason.
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN Rscript setup/install_or_die.r devtools
RUN Rscript setup/install_or_die.r github rexyai/RestRserve
RUN Rscript setup/install_or_die.r languageserver
RUN Rscript setup/install_or_die.r BiocManager

RUN R -e 'BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db"))'

RUN apt-get update && apt-get install --yes libglpk-dev libxt-dev libcairo2-dev
RUN R -e 'devtools::install_github("eddelbuettel/rcppspdlog")'
RUN R -e 'devtools::install_github("kharchenkolab/pagoda2")'
RUN R -e 'devtools::install_github("kharchenkolab/conos@v1.3.1")'

RUN R -e 'install.packages("gprofiler2")'
RUN R -e 'install.packages("RJSONIO")'

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# add app
COPY src ./

# start app
ENTRYPOINT ["bash", "/var/lib/watchfile/entrypoint.sh"]
CMD ["Rscript", "work.r"]

# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

ARG EXPERIMENT_ID
ENV EXPERIMENT_ID $EXPERIMENT_ID

WORKDIR /r
CMD watchmedo auto-restart --directory=./src/ --pattern=* --recursive -- Rscript src/work.r
