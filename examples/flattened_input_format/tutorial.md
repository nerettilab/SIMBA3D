# simba3d Command Line Tools

The simba3d package comes with serveral command line tools to assist in 
managing and looking at results. 

* <'simba3d'> for running simba3d 
* <'simba3d-print'> for viewing specific information from output results
* <'simba3d-disp'> for graphically viewing output results
* <'simba3d-convertion'> for converting from .npz to some other supported format

There is also a python script for each command which can be ran for local 
installs.

* <'python run_simba3d.py'> same as <'simba3d'>
* <'python run_simba3d_print.py'> same as <'simba3d-print'>
* <'python run_simba3d_disp.py'> same as <'simba3d-disp'>
* <'python run_simba3d_convertion.py'> same as <'simba3d-convertion'>

Each command has a command line help option which can be accessed by typing
-h after the command.

This tutorial will outline how to use the command line tools and should be 
followed in succession to work correctly.

## Run simba3d Task list

Start the tutorial in the directory where the tutorial.md file is located.

There are two json task files in the tasks directory which can be used to 
instruct simba3d on what it should run.

The task01.json file is a minimalistic example of a json . In this case, it 
simply tells simba3d where the input data is and where it should output the 
result. All of the settings are set to the default simba3d settings in this 
case.

The task02.json file is also minimalistic, but actually specifies that the
penalties be applied into the optimization.

Run the first task using the simba3d command
> simba3d -r tasks/task01.json

After the program finished running a file will be created with default task name
(which happens to be unnamed_task) appended with a randomly generated string
(a uuid). This can be set manually if desired.

Run the second task using the simba3d command
> simba3d -r tasks/task02.json

Another file will be created in the results directory.

The default names do not convey any information about what experiment was run,
but these files can be loaded into python fairly easily using the numpy module.

## Converting to Other formats

It is recommended the results be saved in the .json since it is more likely to 
load correctly in python. Other output formates can then be converted later 
using the <'simba3d-convertion'> tool.

To convert all the .json files into .mat (MatLab Mat-File format) in the results
directory, type the following command:
> simba3d-convertion results/*.json --ext_out .mat

This command might actually work for generally convering .json to .mat, but this
feature will not be supported to that level of generality. Currently only 
convertion from .json to .mat, .json, and .npz is supported specifically for 
simba3d outputs.

## Printing information from results

The default names of the runs are designed to be unqiue so that it is extreamly 
unlikely that data is overwritten because the same name was picked twice. This
however makes it difficult to distringuish between which output was from which
run.

One way to see the task settings is to use the <'simba3d-print'> command.

To get a list of input keys stored in the file type the following command:
> simba3d-print -i results/output.npz 

To get a specific key value simply append the previous command with a list of
keys you want to print

To get information parameter setting
> simba3d-print -i results/output.npz weight_uniform_spacing weight_smoothing weight_population_prior

We now have a method of determining which result goes to which task.

It can also be used for printing the resulting curve with the following command
>  simba3d-print -i results/output.npz X_evol

One can also pipe it into a file if desired.
>  simba3d-print -i results/output.npz X_evol > result

Unless the user specifally requested not to include the json task file, we can
also get the task information
> simba3d-print -i results/output.npz json_task 

## Displaying graphical summaries of the results

simba3d also come with a graphing tool, if the user has matplotlib installed.

The <'simba3d-disp'> command can be used to display a single result or multiple
results together.

To graph together all the .npz results stored in the results directory use the 
following command:
> simba3d-disp -i results/*.npz

If this tool is working correctly, two figures should appear. Figure 1 shows
penalized energies plotted over the number of iterations. Figure 2 shows the
curves plotted together after removing translation, scale, and rotation.

If no center curve is specified, then all the curves will be centered with 
respect to the first curve in the list. However, if one wishes to control which
one to center to, then one can pass in a results file and center to that one.
> simba3d-disp -c results/output.npz -i results/*.npz

There is also a filter option where the user can filter out parameter settings.
For example if we only want to see results where the weight_smooting parameter 
is zero then we can run the following command:
> simba3d-disp -f weight_smoothing 0 0 -i results/*.npz

This feature is just for getting quick view of results. To get more precise 
controls over the plots, a script should be written specific for your desired
output.
