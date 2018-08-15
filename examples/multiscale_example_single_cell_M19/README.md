# OVERVIEW

This example demonstrates how to create a multi-scale warmstart task on a single
cell dataset. 

The data directory stores the population and the single cell mESC data for
chromosome 19.

# FIRST set up the multiscale warmstart task

The 'create_multiscale_warmstart_tasklist.py' script will set up a multiscale 
warmstart tasklist to estimate the structure. First it will create two sets
of downsampled matrices respectively labeled with the corresponding matrix 
dimensions. Then the script will create a warmstart task list by iteratively
reading the downsampled inputs from the smallest to the largest starting from 
the 70x70 dimesion downsampled matrix all the way up the the 552x552 matrix.
In the first iteration, the initialized curve is randomly generated using
guassian noise. Each subsequence run using the warmstarts will initialize
using the previous iteration's run. To setup the multiscale tasklist, use the 
following command:

> python create_multiscale_warmstart_tasklist.py

# THEN run the tasklist

After running the script, the data directory will contain downsampled matrices
and a tasklist file called 'multiscale_warmstart_tasklist.json' will be made.
This tasklist provides the sequential instructions for SIMBA3D to run the
multiscale warmstarts. To collect experimental results, use the following
command:

> simbda3d -r multiscale_warmstart_tasklist.json

At this point, simba3d will begin collecting the data results and saving the
intermediate runs to the output directory. The intermediate runs are the
solutions to the downsampled verions of the data. The tasklist setup script
'create_multiscale_warmstart_tasklist.py' will randomly create a unique uuid
each time it is run. This uuid is appended to the output name. If you try and 
rerun the experiment without clearing the output, then simba3d will check to 
see if the task has already been collected. It will skip over the task if the 
task in the ouput directoy has the same uuid. This is useful for resuming in the
event of an interuption. The current subtask will be lost, but the intermediate 
tasks that were completed will be kept.

# look at the data

The simba3d module has a display tool for quickly plotting the result collected
from an experiment. This can be done with the following command (the uuid will 
be different because it is randomly generated):

> simba3d-disp -i results/add_a_name_tag_here_[uuid]_70.npz

Two images should pop up. One is the energy evolution and the other is the
estimated curve. The _70 at the end of the file indicates that this result
is created using the 70x70 smallest downsampled matrix. The final result in this
case can be displayed using the following command:

> simba3d-disp -i results/add_a_name_tag_here_[uuid]_552.npz

In this example, there are 7 resamplings of the matrix and the _552 result is
the full scale data.

# alter the tasklist creation script for your own data

The 'create_multiscale_warmstart_tasklist.py' script is documented to act like
a template script for creating multiscale tasklists.

To make one for your own data using this script you may need to edit some parts
of the 'create_multiscale_warmstart_tasklist.py' script. 

Here are the steps:

- Place your data matrices in the 'data' directory
  - You can tell it to look for the data matrices in different directories if 
  you want to, but it can become complicated if you have them all over the 
  place. It is recomended that you put them in one place. You can specify the 
  inputdir (line 25) with a relative or absolut path.
  
- Change the file name specified in the save_matrices function (on lines 32 
and 35) to the name of your matrix files stored in the input directory
  - In this example there are two matrices (a single cell and a population 
  matrix). If you have only one matrix, then comment out the population matrix 
  lines (comment out lines 32, 74, and 109 if you only have one data matrix 
  or if you are not using the population prior term).
  
- optionally, you may wish the change the 'output_filename' (line 76) to 
provide a more intutitive output file name.
  - You might not like to have the uuid in the file name (since it is long). You 
  can edit the output filename so that the UUID string does not appended if you 
  like.
  - Make sure each subtask has a unique name, so that results do not get 
  overwritten.
  
- the 'index_parameters' has some useful preprocessing tools as well. 
  - The offdiagonal entry is set to 1 on this case. In some cases, the 
  offdiagonal entries are unobservable further away from the diagonal, so you 
  may wish to specify that so that simba3d knowns those entries are missing and 
  not observed zeros.
  - The other preprocesssing index tools have a limited amount of testing, but I
  believe they should work.
  
- optionally, you may wish to change the tuning parameters specified in lines
107-109.

- there are other optimization options specified in lines 123-134
  - you can change the maximum number of iterations,
  - whether you want the iterative energies to be displayed while running
  - indicate if you want to store the iterative curves
  - specify the blackbox optimizer method (the defualt is BFGS)
  
- finally you may wish to change the name of the taskfile (line 145)
  - you can provide an absolute or relative file path, but make sure the 
  intermediate directories exist.
  
After editing the script you should be able to setup the multistart tasklist and
run simba3d on the tasklist similiarly as described above.


