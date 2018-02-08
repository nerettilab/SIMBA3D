# About

Structure Inference from Multiscale BAyes in 3 Dimensions (SIMBA3D)

This tool was primarily developed to study the 3D structure Chromatin. It 
includes analytically derived gradients for a variety of penalty functions to
estimate 3D architecture from pairwise interaction matrices.


# Dependencies
## Required Packages
* numpy - needed for scientic programing
* scipy - needed for the optimization tools and other scientic programing

## Recomended Packages
* matplotlib - needed if you want to plot the results in python
* simplejson - recommended over standard json packaged with python


# Installation
After installing the python dependencies, one can install simba3d using 
distutils with the following command

> python setup.py --install

If pip is installed, then it can be installed from this directory using

> pip install .

Pip results in a cleaner install, and it is easier to remove if desired

> pip uninstall simba3d

simba3d can either be used as a python module or as a command line interface.
Using it as a python module has the advantage that a long and tedious sequence 
list of tasks can be generated using python. The command line interface has
the appeal of being easy to run from different directories, so that the user
can create a JSON tasklist in any directory and run it with the `simba3d`
command.

simba3d can also be run locally, but this can limit the file structure 
organization of the project. To do this, copy the 'simba3d' directory (which
contains the '__init__.py') into your project directory. The location of this
directory should be located inside the directory you plan to run python from.
Next create a python script to create the task list, and a python script to
run the simba3d mp_handler.

Examples have been included in the examples directory.

# Command Line Interface

After installing simba3d, the simba3d command should be available.

For help type:
> simba3d -h

If you are using a local install you can also use `python run_simba3d.py` just 
like the `simba3d` command. For example, to get help you can type:
> python run_simba3d.py -h

The help documentation will tell the user that there are two ways to pass in
file parameters. One way is to use the -r option to read in a json file 
containing a tasklist. If a single task is to be ran, then one can pass in input
and output instructions and the parameter settings using the options specified. 
This feature does not currently have all the available parameters and will most
likely be phased out of future versions. It is recommended either individual
files with different tasklists are used or one big file with a large list of 
tasks is used.


## JSON Tasks for simba3d
simba3d requires some information to know what to run. This information can be 
efficiently passed in as a list of tasks stored in a formated json file. 

There are a lot of options available which can make the process of specifying 
the task look complicated at first. There are defaults set for many of the
keys (which are shown in square brackets below).

The individual tasks are dictionaries with the following keys:
* 'taskname' (optional) name of the task [unamed_task]
* 'uuid' (recommended) unique identifier specific to the task [...randomly generated...]
* 'store_task' (optional) store the task in the output [True] 
* 'randomize_initialization' (optional) should we randomize the initialized curve? [False]
* 'seed' (optional) integer between 0 and 4294967295 to fixed seed for randomizing inititialization
* 'check_jacobian' (optional) to check the alytical gradient with numerical gradient [False]
* 'input_file_names' manage file names
    * 'inputdir' parent directory for the input files ['.']
    * "pairwise_contact_matrix"  the pairwise interacation matrix file (path relative to inputdir )
    * "initialized_curve" (optional) file to import initialized curve (path relative to inputdir )
    * "population_contact_matrix" (optional) the pairwise interacation matrix file from poptulation data (path relative to inputdir )
    * "prior_shape_model" (optional) shape prior curve file (path relative to inputdir )
    * 'outputdir' parent directory for the output files ['.']
    * 'output_filename' (optional) the file name to output the result it (.mat or .npz) [uuid.npz]
* 'parameters'
    * 'a' (optional) the a parameter [-3.0]
    * 'b' (optional) not identifiable with unconstrained scale [1.0]
    * 'term_weights' dictionary storing the weights for the penalty terms
        * 'data' (optional) weight for the data term [1.0]
        * 'uniform spacing' (optional)  lambda_1 penalty weight [0.0]
        * 'smoothing' (optional) lambda_2 penalty weight [0.0]
        * 'population prior' (optional) weight for the population matrix prior [0.0]            
* 'options'
    * 'gtol'    (optional) tollerance for the norm of the gradient (stopping criteria) [e-5]
    * 'maxitr'  (optional) set maximum number of iterations [100000]
    * 'display' (optional) display function values at each iteration? [True]
    * 'store'   (optional) store the iterative curves? [False]
    * 'method'  (optional) optimizer option 'BFGS'[default],'AGM','Nelder-Mead','CG','Newton-CG'

You can either write a script to create a list of tasks with the desired 
settings and then pass them into the optimzer, or you can create a JSON file
with a list of tasks that can be read in. JSON is used because it can be 
read and edited manually by a human being and it is flexible to specify many 
settings for a run.

There is also a "data" key which can be used to pass the data matrix and 
initilaized curve as an array (without loading it from a file). This is useful
for running simba3d within a python script. This may be useful for more advanced
users who wish to edit the data matrix (or generate simulated data) without
reading and writing the data to file.

To read a json file with these parameters set, use the following command
> simba3d -r path/to/simba3d_task.json

# How to tell simba3d to run sequentially with warmstarts 

Warmstarts is a technique where the final result from the previous setting is
used to initialized the next run with different settings. This is useful for
parameter tuning because it can lower the total computation time.

The usewarmstarts key must be set to true within a sequentially
ordered task. To create a sequence of task which are to be done in a specific 
order, simply nest them inside an array. In each task inside the array, set
the usewarmstarts to true so that simba3d will know that it should initialize 
the next task with the solution computed from the previous task. 

Several examples can be seen in the examples directory.

# Output File Formats

Currently, can create reports in the following output formats:
* '.npz'  zipped NumPy array file
* '.mat'  MatLab Mat-File file
* '.txt' text file
* '.json' text file

The default format is '.npz' and it is recommended to use this output because
it is easier to load into python and convert to some other format if desired.

A convertion tool has been packaged with simba3d to assist in converting
from .npz to some other desirable format. This convertion tool can be called
from the `simba3d-convertion` command.

For more information on using the convertion tool type
> simba3d-convertion -h

This can also be called on a local install using the command
> python run_simba3d_convertion.py -h

Several examples can be seen in the examples directory.

# Veiwing Results

simba3d has several tools to aid in summarizing the results collected.

The `simba3d-print` command can be used to get print out specific key data
values stored in the results file. 

If the matplotlib module is installed, then simba3d also has a graphical
summary command which can be called by `simba3d-disp`

This command will display the final curve and the energy of selected results.

Several examples can be seen in the examples directory.
