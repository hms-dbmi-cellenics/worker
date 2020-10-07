#!/bin/bash

# Bash "strict mode", to help catch problems and bugs in the shell
# script. Every bash script you write should include this. See
# http://redsymbol.net/articles/unofficial-bash-strict-mode/ for
# details.
set -euo pipefail

# Install required dependencies for Rserve.
export R_CRAN_PKGS="Rcpp R6 uuid checkmate mime jsonlite"
xargs Rscript setup/install_or_die.r $R_CRAN_PKGS

# Install Rserve.
Rscript -e "install.packages('Rserve', repos = 'http://www.rforge.net/')"

# Install RestRserve.
Rscript -e "remotes::install_github('rexyai/RestRserve')"