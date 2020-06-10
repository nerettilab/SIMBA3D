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

> simbda3d -i multiscale_warmstart_tasklist.json

At this point, simba3d will begin collecting the data results and saving the
intermediate runs to the output directory. The intermediate runs are the
solutions to the downsampled verions of the data. The tasklist setup script
'create_multiscale_warmstart_tasklist.py' will randomly create a unique uuid
each time it is run. This uuid is appended to the output name. If you try and 
rerun the experiment without clearing the output, then simba3d will check to 
see if the task has already been collected. It will skip over the task if the 
task in the ouput directoy has the same uuid. This is useful for resuming in the
event of an interuption. The tasks that were completed will be skipped when you
resume running.

