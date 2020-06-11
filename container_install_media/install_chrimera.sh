#!/bin/bash

# install wget
apt-get update \
  && apt-get install -y wget
  
# Download chimera
wget https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64/chimera-1.14-linux_x86_64.bin
chmod +x chimera-1.14-linux_x86_64.bin

# run chimera install command
./chimera-1.14-linux_x86_64.bin