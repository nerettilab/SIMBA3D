# Make a pdb file with the output

A pdb file is an input for UCSF Chimera. It is an efficient visualizing tool 
written in C (or C++?) which is used for visualizing molecules.

You can create a .pdb file by placing setting the 'make_pdb_file' to true
in the json task. This will create a similiarly named pdb file.

In the tasks directory there are two tasks. One has 'make_pdb_file' set to 
false, and the other set to true. Run the two tasks to see what happens with 
the following command.

> simba3d -i ./tasks/*

This should create three files.

If you ran a result without this flag decide you want pdb file, you do not need
to rerun the experiment. You can use the following command to create pdb files
for all the results.
> simba3d-convertion --ext_out .pdb results/*.json
