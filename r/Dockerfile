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

# install pip and Radian for interactive R shell
RUN pip install -U jedi radian

# install required R packages for the worker.
COPY requirements_r.txt .
RUN cat requirements_r.txt | tr '\n' '\0' | xargs -0 -n1 xargs Rscript setup/install_or_die.r && rm -r setup

# add app
COPY src ./

# start app
CMD ["Rscript", "work.r"]