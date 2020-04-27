# pull official base image
FROM python:3.7-buster

# set working directory
WORKDIR /src

# install app dependencies
COPY requirements.txt ./
RUN pip3 install -r requirements.txt

# add app
COPY src ./

# start app
CMD ["python3", "-u", "work.py"]