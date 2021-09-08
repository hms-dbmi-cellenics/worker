# pull official base image
FROM python:3.8-buster AS builder

# set working directory
WORKDIR /src

# install app dependencies
COPY requirements.txt ./
RUN pip3 install -r requirements.txt


# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# add app
WORKDIR /src
COPY src ./

# start app
ENTRYPOINT ["bash", "/var/lib/watchfile/entrypoint.sh"]
CMD ["python", "-m", "worker"]

# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

RUN pip3 install -U watchdog[watchmedo]

WORKDIR /python/src
CMD watchmedo auto-restart --directory=. --pattern=*.py --recursive -- python -m worker
