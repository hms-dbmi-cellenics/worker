# pull official base image
FROM rexyai/restrserve:0.3.0

# set working directory
WORKDIR /src

# install app dependencies
COPY install_packages_or_die.R ./
COPY requirements.txt ./
RUN cat requirements.txt | xargs Rscript install_packages_or_die.R 

# add app
COPY src ./

# start app
CMD ["Rscript", "work.r"]