# OVERVIEW

Creatings long json tasklists for simba3d by hand is tedious. One way
to streamline this process is to use a python script to create the
json file.

Another advantage over manually editing a json file is that a python script 
allows comments so that a template for creating different types of tasks can
be made and more easily modified. 

This also makes it easier to create more complicated file structures so that
thousands of results do not flood the working directory.

With this approach, first a python script is run to create a json tasklist, and
then the <'simba3d'> command is run on that tasklist.

It is recomended that simplejson is installed since it is more likely to be
updated than the standard json package.

# simple but long tasklist example

A python script called create_simple_but_long_tasklist.py has been created
which illustrates how to create an arbitrarily long tasklist using python. In
this example, a search over various settings of uniform_spacing penalty is
specified. More specifically, lambda_2 = lambda_3 = 0.5 but lambda_1 varied
from 0 ,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0.

To create this tasklist run to following command:
> python3.6 create_simple_but_long_tasklist.py

After running the script, a file called simple_but_long_tasklist.txt will be
made. Once this tasklist is created, one can run it using the following command:
> simba3d -i simple_but_long_tasklist.json

# warmstarts tasklist example

A python script called create_warmstarts_tasklist.py has been created
which illustrates how to create a warmstarts tasklist.This is accomplished
by setting the 'usewarmstarts' parameter to True and nesting the tasklist inside
another list. As before a search over various settings of uniform_spacing 
penalty is specified. More specifically, lambda_2 = lambda_3 = 0.5 but lambda_1 
varied from 0 ,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0. simba3d will 
run the list of task sequentially, and if 'usewarmstarts' is set to true it 
will use the solution from the previous run to initialize the next.

To create this tasklist run to following command:
> python3.6 create_warmstarts_tasklist.py

After running the script, a file called warmstart_tasklist.txt will be
made. Once this tasklist is created, one can run it using the following command:
> simba3d -r warmstart_tasklist.json

# compare the two

You can create a latex report with the following command
> simba3d-result-disp -p with_and_without_warmstarts -i results/* 

This should create an ordered report with the two set of runs. One will have 
the warmstarts run which uses the previous run to initialize the next, and 
the other will have a random initialization for each run.

# run simba3d within python script example

simba3d can also be called within the script itself by loading the mp_worker.

An example of this can be seen in the run_simba3d_within_python_script.py 
script. Run the script with the following command:
> python3.6 run_simba3d_within_python_script.py

This time, the simba3d command does not need to be called because the script
will run the tasklist.
