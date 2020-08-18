In this example, you will find multiple input formats that simba3d can read in
order to specify the pairwise interaction matrix.

# Introduction

There are two main input format types for data matricies to input into SIMBA3D.

There is a dense matrix representation which is inputed as an array of row 
arrays. This is not ideal since the data is often sparse and it is symmetric
for inter cromosomal data.

There is also a sparse representation which is inputed as three separate 
arrays which store the row index, column index, and value. The sparse input 
format is very useful for single cell, and it is closer to raw data formats.

We support both types of inputs for convience.

# json, mat, npy, and csv

There are almost as many formats for passing data as there are data sets. Some
data formats are binaries which require knowing how to read the header of the
file in order to read the data. Matlab's .mat file is a binary. Python's .npy
files are binaries as well. On the other hands, json is a plain text format 
that can be read using just about any code language. A csv (comma separate 
value) consists of data values separated by one or more delimiters (usually
a comma and an end line character).

We suggest json or csv as a default format to passing in data into SIMBA3D, so that 
users are not restricted to use python or matlab for data manipulations. If you 
choose to use mat or npy files as inputs, then anyone who wants to use
your data will be required to use python or matlab if they want to reprocess it.
If you use json or csv, then most people will be able to read it without much
assistance.

The csv format is probably more familiar to biologist and statisticians. It is 
less flexible than json though.

The sections below go into further details.

# Considerations for Future Format Support

npy and npz formats may not be supported in future versions because the binaries
are not compatible between python 2 and python 3. This is a huge problem if you
have multiple systems sharing a workload, or if you decide to upgrade to python
3 since python 2 is no longer supported. For these reasons, the native python
data formats are no longer recommended for storing chromosome reconstruction
data with SIMBA3D. There are no plans to actively disable them, but if they
happen to break, it will be a low priority.

If you choose to use the matlab binary, then please note that matlab is
proprietary software, and they could choose to change the binary format in
future releases of matlab so that they cannot be read by the python scipy.io
package. As a result, SIMBA3D might not be able to read such a file.

# Sparse Matrix Input Formats for SIMABA3D

The "sparse_pairwise_contact_matrix" and "sparse_population_contact_matrix" 
input file can read various formats.

The sparse representation consists of three arrays called "row_index", 
"column_index", and "count". The three arrays must have the same length. We also 
need the row_dimension and the column_dimension specified as integers. For 
intochromosomal data, we only need the upper triangular (or lower triangular)
part.

For those of you who are familiar with matlab, you can just create a mat file
containing these 5 items and SIMBA3D will be able to read that as a sparse
matrix input.

The default input for SIMBA3D is a json file containing these five keys. This
is chosen for the default because it is human readable and does not center
around a proprietary binary format.

Both the mat and json formats have an advantage that additional information
can be easily tagged along with the data.

# Dense Matrix Input Formats for SIMABA3D

The "pairwise_contact_matrix" and "population_contact_matrix" input file can 
read various formats. The data can be stored as either a npy, npz, matlab, or 
json file.

For .mat, simply create a matrix and name the variable either 
"pairwise_contact_matrix" or "population_contact_matrix".

For .json, simply create a nested array (an array containing an ordered list
of arrays corresponding the the rows of a matrix) and name the variable either 
"pairwise_contact_matrix" or "population_contact_matrix".

For .csv, simply create a comma separated data matrix with column entries
separated by commas and row entries separate by end line characters.

# Experiment

Within the data_inputs directory there are various input formats which can
be read by SIMBA3D. In the tasks sub-directory there are corresponding tasks
which load each corresponding format type. You can run them with the following
command. 

> simba3d -i tasks/*

You should call this command from the directory which contains this readme file
since the tasks use relative paths.

This will begin to reconstruct things using various input examples.

You should get 6 identical results with different computation times. Use the 
following command to generate a report to compare the results. 
> simba3d-result-disp -p formats -i ./*.json


# Conclusion

There are many option and considerations for how you choose to interface with
SIMAB3D.
