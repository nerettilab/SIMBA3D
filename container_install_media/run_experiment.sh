#!/bin/bash

cd experiments/$1
echo "Running Tasks"
simba3d -i tasks/*.json
