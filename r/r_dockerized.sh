#!/bin/sh -
docker image build -t worker-r .
docker run worker-r --name worker-r-dev
docker exec -it worker-r-dev R "$@"