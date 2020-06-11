#!/bin/bash

# install git
#apt-get update \
#  && apt-get install -y git

# clone the repository
#git clone https://github.com/nerettilab/SIMBA3D.git

tar -zxvf SIMBA3D-3.0.0.tar.gz

ls -lah
# move into the cloned repository
cd SIMBA3D-3.0.0
ls -lah
# use pip3 to install it.
pip3 install .

# create mounts
mkdir -p /experiment/data
mkdir -p /experiment/tasks
mkdir -p /experiment/results
mkdir -p /experiment/reports
