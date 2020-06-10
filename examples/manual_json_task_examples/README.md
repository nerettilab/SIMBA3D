# OVERVIEW
In this folder there are various json files which illustrate how to pass 
instructions into simba3d for a variety of tasks. This tutorial is for users
who wish to edit the json file manually to instruct simba3d to run a variety
of tasks.

# minimalist example

This example contains the bare minimum required to get simba3d to run something.

The minimal.json file contains a task list with a single task that tells
simba3d which pairwise_contact_matrix to load.

To run the task list use the following command
> simba3d -i minimal.json

Each time you run the command above, simba3d will create a new json file that 
has a name that looks like
'''
unnamed_task_[long random string of letters and numbers].json
'''
This happens for several reasons. 
1) The 'file_names_output_filename' which allows the user to specify name of the 
output file was not specified in the task json file. When this happens simba3d defaults to generating a name by appending the taskname and uuid.
	- You only set the 'file_names_output_filename' if you do not like the way simba3d automatically names the output file using the taskname and uuid.
2) The 'taskname' was not specified so simba3d sets the taskname to 'unnamed task' by default.
3) The 'uuid' was not specified, so simba3d defaults to randomly generating one for that task.

simba3d will check if the outputfile specified by task already exists and if it has the same uuid it will skip over that task because it is presumes to have already collected that task. 
This is useful for resuming long task lists that have been interrupted. 
Since 'file_names_output_filename', 'taskname', and 'uuid' are not specified, the outputfile name and uuid are randomly generated and simba3d has no way
of telling that the task has already been ran.
This means that it will almost certainly not find a completed task with the same uuid, and it will rerun the same exact task each time you run this command.


# minimal recommended example

While simba3d has defaults for most parameters, it is not recommended to 
relay on them because future versions may change these defaults and it obscures
what was actually run. 
Also, it makes it easier to modify the task json if the keys are already there.
It is recommended at a minimum that the user specifies 
the keys stored in the minimal_recommended.json file. 

It is recommended that the user provide a short descriptive name for the task 
and provide a uuid if it is undesirable to rerun this task more than once. An 
optional more detailed description is also nice to add since it makes the
task more readable.

It is also recommended that the user provide an initialized curve, otherwise
simba3d will initialize with the default initialization (which may change in
future versions).

To run the task list use the following command
> simba3d -i minimal_recommended.json

Notice that the outputfile will be defined using the taskname and the uuid. If
the task is rerun, it will look for a task with the same name and check the 
uuid. It will skip the experiment if it has already been collected.

# multiple independent tasks example

This task is a remnant for a feature which is not supported and may become 
deprecated. I am talking about the multi-processing feature which runs tasks 
in parallel. The feature may still work, but I will not be supporting it 
because I have found it is too difficult to debug and to get it working generally across 
platforms.

The first layer of the json task list is reserved for tasks that are to be 
completed independently of other tasks. This is most often the case unless
the user wants to control the order that the tasks are completed in when using
the multiprocessing feature or if sequential warmstarts being done to reduce
computation time.

To run the task list use the following command
> simba3d -i multiple_independent_tasks.json

# sequential warmstarts example

If you compare the multiple_independent_tasks.json and sequential_warmstarts.json,
you will see that the main difference is that the two tasks are lumped into one 
bracket. 
This tells simba3d, the task are related and should be completed 
sequentially in the order provided. 
You should also note the 'usewarmstarts' key is set to true.

The main reason for allowing sequential tasks is to place a mechanism for 
performing warmstarts. Warmstarts is a technique typically used to speed up
parameter tuning searches by using the solution from the previous setting
to initialize the solution for the next parameter setting. If the solution
is expected to change slightly for small changes in the parameter, this
can reduce the time needed to do a grid search considerably. There is also a
risk that the solution space may not be fully searched by using this approach.

As it is, to tell simba3d to use warmstarts inside a sequential task list,
simply set the 'usewarmstarts' key to true inside the first nested sequential 
task list.

To run the task list use the following command
> simba3d -i sequential_warmstarts.json
