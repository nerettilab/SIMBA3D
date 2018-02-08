# OVERVIEW
In this folder there are various json files which illustrate how to pass 
instructions into simba3d for a variety of tasks. This tutorial is for users
who wish to edit the json file manually to instruct simba3d to run a variety
of tasks.

# minimalist example

This example contains the bare minimum required to get simba3d to run something.

The minimal.json file contains a task list with a single task that points tells
simba3d which pairwise_contact_matrix to load.

To run the tasklist use the following command
> simba3d -r minimal.json

Each time you run the command above, simba3d will create a random uuid for the 
task. The uuid is a unique string which is specific to a task. SIMBA3D will 
search the output directory for task and if it finds one with the same uuid it
will skip over it because it is presumed to be already ran. This is useful for 
resuming long tasklists. Since no uuid is provided, simba3d will almost 
certainly not find a completed task with the same uuid so it will rerun the
same exact task.

Currently, the curve is initialized randomly, but future versions may use some
other initilization by default.

# minimal recommended example

While simba3d has defaults for most parameters, it is not recommended to 
relay on them because future versions may change these defaults and it obsures
what was actually run. It is recommended at a minimum that the user specifies 
the keys stored in the minimal_recommended.json file.

The minimal_recommended.json file contains a list with a single task. Even if
the default options and parameters are used for this task it is still 
recommended to specify them to more easily edit the task manually.

It is recomended that the user provide a short descriptive name for the task 
and provide a uuid if it is undesirable to rerun this task more than once. An 
optional more detailed description is  also nice to add since it makes the
task more readable.

It is also recommended that the user provide an initialized curve, otherwise
simba3d will initialize with the default initialization (which may change in
future versions).

To run the tasklist use the following command
> simba3d -r minimal_recomended.json

Notice that the outputfile will be defined using the taskname and the uuid. It
the task is rerun, but a task with the same uuid is found, simba3d will skip 
over the task so that it will not be rerun.

# multiple inpendent tasks example

The first layer of the json tasklist is reserved for tasks that are to be 
completed independently of other tasks. This is most often the case unless
the user wants to control the order that the tasks are completed in when using
the multiprocessing feature or if sequential warmstarts being done to reduce
computation time.

To run the tasklist use the following command
> simba3d -r multiple_independent_tasks.json

simba3d is setup to run embarassingly parallel tasks using the multiprocessing
python module. To use this feature use the -c option on the command line. For
example, if you want to run the task in parallel two at a time use the following
command.
> simba3d -c 2 -r multiple_independent_tasks.json

This is essentially equivalent to opening two separate instances of python and
running the first and second task at the same time. When a task is completed, 
simba3d will automatically look for the next task in the list.

# sequential tasks example

By nesting a list of tasks within a list, that list then becomes one task with
multiple subtacts which are to be completed. In this example the tasks are
nested inside of a list which tells simba3d that these are to be completed
sequentially in the order provided.

To run the tasklist use the following command
> simba3d -r sequential_tasks.json

Notice that if you tell simba3d to run the task 2 at a time, that there is no
boost in speed. This is because there is only one independent task in the list
which contains two tasks that are to be completed sequentially.

# sequential warmstarts example

The main reason for allowing sequential tasks is to place a mechanism for 
performing warmstarts. Warmstarts is a technique typically used to speed up
parameter tuning searches by using the solution from the previous setting
to initiailize the solution for the next parameter setting. If the solution
is expected to change slightly for small changes in the parameter, this
can reduce the time needed to do a grid search considerably. There is also a
risk that the solution space may not be fully searched by using this approach.

As it is, to tell simba3d to use warmstarts inside a sequential tasklist,
simply set the 'usewarmstarts' key to true inside the first nested sequetial 
tasklist.

To run the tasklist use the following command
> simba3d -r sequential_warmstarts.json

# supported output filename formates example

There are four file formats that the data can be outputed to (.npz, .mat, .txt,
and .json). To specify the file format for the output file, set the 
"output_filename" key nested within the "file_names" dictionary. Append the
desired extention to the file to ouput in a different format. The default file
format is the zipped NumPy array file '.npz' which is selected if no extention
is provided.

To run the tasklist use the following command
> simba3d -r supported_output_filename_formates.json

This will run separate experiments for each of the output formats.

Currently there is no way to output multiple formats in one task. Also, simba3d
is currently only set up to read npz and mat files from outputs. If other 
formats are used, simba3d will not be able to find tasks that have already been
run to resume an incomplete tasklist. This may eventually be addressed in 
future verions. It is recommended that the .npz or .mat format 
is used and then converted to desired format using a python script. Some 
examples have been included with the source module.
