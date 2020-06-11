#!/bin/bash

cd experiments/$1
echo "Generating Report"
simba3d-result-disp -o reports -p sum -i results/*.json
