# Guidance of Structuring the Experiments Directory

1) Put you raw data into the raw_data directory.
2) Edit the your_raw_data_proccessing_script.sh to read the raw data and convert it into something that simba3d can read and strore it in the data directory.
3) For each experiment, create a directory that  has reports, results, and tasks.
	- Set "file_names_inputdir"  to "../data/whatever-subfolder" in your tasks files
	
# Running in a docker container.

The following command will run the simba3d image and execute an interactive bash terminal
with your experiments directory mounted
> docker run --rm -it -v /absolute/path/to/experiments:/experiments simba3d bash

Use the cd command to move into the directory
> cd experiments/experiment01

Then you can use the simba3d commands like usual.

To run all the tasks
> simba3d -i tasks/*.json

To create a report of the tasks
> simba3d-result-disp -o reports -p summary -i results/*.json 

