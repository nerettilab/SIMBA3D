From python

# copy the script which will generally run the experiment
COPY ./container_install_media/*.sh ./

# copy the setup script
COPY ./dist/SIMBA3D-*.tar.gz ./

# give execute permissions
RUN chmod +x *.sh

# remove any windows newline characters
RUN sed -i 's/\r//g' ./install_simab3d.sh

# run the install script
RUN ./install_simab3d.sh
