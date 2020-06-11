#!/bin/bash

rm -rf ./dist/*
python3 setup.py sdist
# if you are running simba3d in docker, use this command
sudo docker build . -t simba3d
