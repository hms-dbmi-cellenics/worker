#!/bin/bash

# Bash "strict mode", to help catch problems and bugs in the shell
# script. Every bash script you write should include this. See
# http://redsymbol.net/articles/unofficial-bash-strict-mode/ for
# details.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

# Install required dependencies for Rserve.
apt-get update
apt-get -y upgrade
apt-get -y install --no-install-recommends r-cran-rcpp r-cran-uuid r-cran-r6 r-cran-checkmate r-cran-mime r-cran-jsonlite

# Install Rserve.
Rscript -e "install.packages('Rserve', repos = 'http://www.rforge.net/')"

# Install RestRserve.
Rscript -e "remotes::install_github('rexyai/RestRserve')"