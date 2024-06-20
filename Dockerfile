FROM ubuntu:22.04

RUN apt-get update && apt-get -y upgrade

RUN apt-get install -y python3-pip python3.10-venv wget pkg-config libhdf5-dev

RUN mkdir /home/user
COPY requirements.txt /home/user/requirements.txt

RUN python3 -m venv /home/user/venv

# ARG を使用して NUPACK_PK を受け取る
ARG NUPACK_PK
COPY "${NUPACK_PK}" "/home/user/${NUPACK_PK}"
RUN /home/user/venv/bin/pip install "/home/user/${NUPACK_PK}"

RUN /home/user/venv/bin/pip3 install -r /home/user/requirements.txt

RUN apt-get install -y cmake make 

# CMD ["/bin/bash"]
# ENTRYPOINT ["/home/user/venv/bin/python3"]

