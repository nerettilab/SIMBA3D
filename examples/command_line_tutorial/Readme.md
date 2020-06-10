# simba3d Command Line Tools

The simba3d package comes with serveral command line tools to assist in 
managing and looking at results. 

* `simba3d` for running simba3d 
* `simba3d-print` for viewing specific information from output results
* `simba3d-results-disp` for graphically viewing output results
* `simba3d-convertion` for converting from .npz to some other supported format

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

Another file will be created in the results directory. These files can be 
loaded into python fairly easily using the numpy module.

## Converting to Other formats

It is recommended the results be saved in the .json format. Other output 
formates can then be converted later using the `simba3d-convertion` tool.

To convert all the .npz files into .mat (MatLab Mat-File format) in the results
directory, type the following command:
> simba3d-convertion results/*.npz --ext_out .mat

This command might actually work for generally convering .npz to .mat, but this
feature will not be supported to that level of generality. Currently only 
convertion from .npz to .mat, .json, and .txt is supported specifically for 
simba3d outputs.

It can also be used to output .pdb file to be read by UCSF Chimera, but this 
outputed format is only for display purposes (it does not store any infromation 
about the run).

## Printing information from results

The default names of the runs are designed to be unqiue so that it is extreamly 
unlikely that data is overwritten because the same name was picked twice. This
however makes it difficult to distringuish between which output was from which
run.

One way to see the task settings is to use the `simba3d-print` command.

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


## Displaying graphical summaries of the results

simba3d also come with a graphing tool, if the user has matplotlib installed.

The `simba3d-disp` command can be used to display a single result or multiple
results together.

To graph all the results stored in the results directory use the 
following command:
> simba3d-disp -i results/*

