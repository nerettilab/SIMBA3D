# OVERVIEW

This example demonstrates how to create a multi-scale warmstart task on a 
simulated data set.

This tutorial is similiar to the 'multiscale_example_single_cell_M19' example
but there is the added step of creating a simulated data matrix from a ground
truth curve stored in the 'data/ground_truth_curves/double_spiral_curve.npy'.

It is useful to use a simulated data set because in real situations the 
groundtruth shape of the 3D structure is unknown. This allows the investigator
to judge how robust the algorithm is to noise and other data artifacts.

# FIRST create a simulated data matrix

> python3.6 create_simulated_data.py

The 'create_simulated_data.py' script will read the 3D curve stored in 
'data/ground_truth_curves/double_spiral_curve.npy' and will resample it to the
desired number of bins and create a simulated pairwise interaction matrix using
independent Poission random variables. The means are computed using the 
Varoquaux model. 

The script will create images of the ground truth curve, and the simulated
data matrix. The data matrix plot label non-zero simply shows which entries are
greater than zero. This helps gauge the amound of sparcity in the simulated 
data.

# NEXT set up the multiscale warmstart task

The 'create_multiscale_warmstart_tasklist.py' script will set up a multiscale 
warmstart tasklist to estimate the structure. First it will create a set
of downsampled matrices respectively labeled with the corresponding matrix 
dimensions. Then the script will create a warmstart task list by iteratively
reading the downsampled inputs from the smallest to the largest.

In the first iteration, the initialized curve is randomly generated using
guassian noise. Each subsequence run using the warmstarts will initialize
using the previous iteration's run. To setup the multiscale tasklist, use the 
following command:

> python3.6 create_multiscale_warmstart_tasklist.py

# THEN run the tasklist

After running the script, the data directory will contain downsampled matrices
and a tasklist file called 'multiscale_warmstart_tasklist.json' will be made.
This tasklist provides the sequential instructions for SIMBA3D to run the
multiscale warmstarts. To collect experimental results, use the following
command:

> simbda3d -i simulated_data_multiscale_warmstart_tasklist.json

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

# viewing the results

You can create a latex report with the following command:

> simba3d-disp -c data/ground_truth_curves/double_spiral_curve.npy  -p sim -i results/*.json 

