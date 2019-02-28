# Overview

For larger matrices, the BFGS method may take up too much memory and may 
result in the use of swap memory. This means that your system might start 
reading and writing to file, which is very slow. This happens with large high 
dimensional optimizations.

Here are two ways to mitigate this problem:
(1) Increase the amount of memory
(2) Try the limited memory version of BFGS

The limited memory version of BFGS uses less stored iterations to estimate the
hessian.

If you are estimating a very large curve, then I recommend using the limited 
memory version. Otherwise, you may find that it takes hours to complete a 
single iteration.

To use the limited memory BFGS, simply set the method to 'L-BFGS-B'. 

# running tutorial

First create an example data matrix using the following command:

> python create_simulated_data.py

This will make an huge data matrix. 

There are two tasks in the tasks folder. One uses 'BFGS' and the other uses
L-BFGS-B'. Run the one with L-BFGS-B with the following command:

> simba3d -i tasks/000001_using_L-BFGS-B_tasklist.json

It will slowly iterate, but will eventually complete the large task within a
reasonable time frame.

Just for fun, you can try using the standard BFGS with the following command:

> simba3d -i tasks/000001_using_BFGS_tasklist.json

It will probably iterate much slower because it uses up more memory and is 
more likely to use swap memory (which reads and writes to disk).

