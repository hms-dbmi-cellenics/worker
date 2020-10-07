#!/bin/bash

# Bash "strict mode", to help catch problems and bugs in the shell
# script. Every bash script you write should include this. See
# http://redsymbol.net/articles/unofficial-bash-strict-mode/ for
# details.
set -euo pipefail

while IFS= read -r line || [[ -n "$line" ]]; do
    Rscript setup/install_or_die.r $line
done < requirements_r.txt