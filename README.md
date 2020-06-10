# About

Structure Inference from Multiscale BAyes in 3 Dimensions (SIMBA3D)

This tool was primarily developed to study the 3D structure Chromatin. It
include analytically derived gradients for a variety of penalty functions to
estimate 3D architecture from pairwise interaction matrices.

# What is New

I am deprecating multi-processing part so that only the solo run
will work.
I have found the multi-processing messages too difficult to debug,
and I just cannot support it.
Sorry, if you use it.
You can still use it in as a python script, but I will not be updating the 'simba3d_taskrun.py' file
with new features.

I have improved the input json to something easier for a human being to input.
See the section below on Simplified JSON format.
I do not plan to purposely remove the old format, but I will not be documenting it anymore.
I think most people will agree that it is an improvement.

The results generated can now be converted to a pdb file which can be viewed in
UCSF Chrimera.

There are now tools for generating automated reports.

The examples have been updated.

# What is coming up

I plan to implement a cythonized version of the code to increase the computation
speed. I will be putting that in a separate branch and I plan to maintain a pure
python version of everything going forward. What this means is that there will
be c/c++ code in the source which will need to be built for your specific
system in order to run if you want to use the cythonized branch. This make the
installation a little more involved, which I why I want to put them in separate
branches.

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

If you do not have root access you can do an install using the user scheme

pip install . --user

This is probably safer to use in general. Installing using root privelges could
potentially create os issues in some cases. It is unlikely that the os uses
numpy or scipy or any other packages simba3d may want updated.

simba3d can either be used as a python module or as a command line interface.
Using it as a python module has the advantage that a long and tedious sequence
list of tasks can be generated using python. The command line interface has
the appeal of being easy to run from different directories, so that the user
can create a JSON tasklist in any directory and run it with the `<simba3d>`
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

If you are using a local install you can also use <'python run_simba3d.py'> just
like the <'simba3d'> command. For example, to get help you can type:
> python run_simba3d_solo.py -h

## Simplified JSON format

All the previous formats should still work, but the json format from here on out has been flattened.
This should be easier for the user because the nesting brackets can be complicated to keep track of, especially when cutting and pasting examples.

The individual tasks are dictionaries with the following keys:
* 'taskname' (optional) name of the task [uuid]
* 'uuid' (recommended) unique identifier specific to the task [randomly generated]
* 'store_task' (optional) store the task in the output [True]
* 'randomize_initialization' (optional) should we randomize the initialized curve? [False]
* 'seed' (optional) integer between 0 and 4294967295 to fixed seed for randomizing inititialization
* 'check_jacobian' (optional) to check the alytical gradient with numerical gradient [False]
* 'file_names_inputdir' parent directory for the input files ['.']
* 'file_names_pairwise_contact_matrix"  the pairwise interacation matrix file (path relative to inputdir )
* 'file_names_sparse_pairwise_contact_matrix"  the pairwise interacation matrix file in a sparse format (path relative to inputdir )
* 'file_names_initialized_curve" (optional) file to import initialized curve (path relative to inputdir )
* 'file_names_population_contact_matrix" (optional) the pairwise interacation matrix file from poptulation data (path relative to inputdir )
* 'file_names_sparse_population_contact_matrix" (optional) the pairwise interacation matrix file from population data in a sparse format (path relative to inputdir )
* 'file_names_prior_shape_model" (optional) shape_prior curve file (path relative to inputdir )
* 'file_names_outputdir' parent directory for the output files ['.']
* 'file_names_output_filename' (optional) the file name to output the result it (.json, .mat, or .npz) [taskname_uuid.json]
* 'parameters_a' (optional) the a parameter [-3.0]
* 'parameters_'b' (optional) not identifiable with unconstrained scale [1.0]
* 'parameters_term_weights' dictionary storing the weights for the penalty terms
* 'parameters_term_weights_data' (optional) weight for the data term [1.0]
* 'parameters_term_weights_uniform_spacing' (optional)  lambda_1 penalty weight [0.0]
* 'parameters_term_weights_smoothing' (optional) lambda_2 penalty weight [0.0]
* 'parameters_term_weights_population_prior' (optional) weight for the population matrix prior. You also need to input the population matrix for this [0.0]        
* 'parameters_term_weights_pairwise_repulsion' (optional) weight for the term which keeps points from getting to close [0.0]    
* 'options_gtol'    (optional) tollerance for the norm of the gradient (stopping criteria) [e-5]
* 'options_maxitr'  (optional) set maximum number of iterations [100000]
* 'options_display' (optional) display function values at each iteration? [True]
* 'options_store'   (optional) store the iterative curves? [False]
* 'options_method'  (optional) optimizer option 'L-BFGS_B'[default],'BFGS','AGM','Nelder-Mead','CG','Newton-CG'
* 'index_parameters_missing_rowsum_threshold'   (optional) specifies a threshold matrix row sum to treat an entry as missing data (missing nodes are ignored in the optimization)
* 'index_parameters_index_missing'              (optional) specifies which entries are treated as missing (missing nodes are ignored in the optimization)
* 'index_parameters_off_diagonal'               (optional) offset off of the diagonal entries which is treated as missing (and ignored)
* 'index_parameters_pairwise_combinations'      (optional) optionally specify specific pairwise combination to use
* 'make_pdb_file'   (optional) [false] boolean which states if you want to output a pdb (to view with UCSF Chrimera) file along with the outputed results

Examples can be seen in examples/flattened_input_format

# How to tell simba3d to run sequentially with warmstarts

Warmstarts is a technique where the final result from the previous setting is
used to initialized the next run with different settings. This is useful for
parameter tuning because it can lower the total computation time. It can also
limit space of curves you can converge to, which may not always be ideal since
the previous setting may trap the curve in a local solution.

In any case, the usewarmstarts key must be set to true within a sequentially
ordered task. To create a sequence of task which are to be done in a specific
order, simply nest them inside an array. In each task inside the array, set
the usewarmstarts to true so that simba3d will know that it should initialize
the next task with the solution computed from the previous task.

Several examples can be seen in the examples directory.

# Output File Formats

The default format is no longer '.npz'. The python2 and python3 binaries are not compatable.
The default output is now '.json'. You can still output to a .npz or .mat file.
