# pull official base image
FROM r-base:4.0.2

# set working directory
WORKDIR /src

# These are the default packages that we must have installed
# for Rserve and RestRserve. They will not change so we want to
# install them at the beginning of the image.
COPY setup ./setup

# install required debian packages to make R packages work.
COPY requirements_debian.txt .
RUN cat requirements_debian.txt | xargs ./setup/install_debian_packages.sh

# install Radian for interactive R shell
RUN pip install -U jedi radian

# install required R packages for the worker.
# these are sectioned as seen here so we can cache compiled code
# NOTE: add new packages to the END of this list. this makes sure
# previously installed packages don't get recompiled for no reason.
RUN Rscript setup/install_or_die.r devtools
RUN Rscript setup/install_or_die.r github rexyai/RestRserve
RUN Rscript setup/install_or_die.r languageserver

# add app
COPY src ./

# start app
CMD ["Rscript", "work.r"]