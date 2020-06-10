# Overview

You can initialize your curves for SIMBA3D and write them to file as either a
csv or a json format. This gives the experimentor more control and the ability
to reproduce results without using a fixed seed.

To aid in creating initializations in the proper format, an initialization 
script has been provided which will create randomly generated initializations 
which can be loaded by SIMBA3D.

# Review the help file

Use the following command to see the help file: 

> simba3d-initialize.py -h

# Create an initialized curve with 100 points

> simba3d-initialize -n 100 -o inputs/initialized_curves/initialization00.csv

# Create several initialized curves

> simba3d-initialize -n 100 -o inputs/initialized_curves/initialization01.csv initialized_curves/initialization02.csv initialized_curves/initialization03.csv

# Make hundreds of intializations using the command

> simba3d-initializeze -n 100 -o inputs/initialized_curves/initialization{0001...1000}.csv 

# initialize two tasks with precisely the same curve

An example task has been provided in the tasks folder. To run the task use the
following command.

> simba3d -i tasks/*

Since there is no output or uuid provided for the task, it will create a 
randomly generated uuid for the output file, and it will not be able to
determine if the result has been ran. This means that if you rerun the
same command it will not think the task has been completed and it will rerun. 

Go ahead and retype the same command.

> simba3d -i tasks/*

After the exiperiment finished running you will find two results with identical
outputs (except for the uuid and file name).

Create a report with the following command:
> simba3d-result-disp -p same -i *.json

You will see that everything except the name and computation time are the same.
