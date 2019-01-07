# Sparse Matrix Format

A sparse matrix data file can be created using a json format. 
The non-empty entries of the matrix are specified using three arrays
corresponding to the row index, column index, and count. The three arrays
must have the same number of elements in them which corresponds to the
number of observable interaction pairs. You must also specify the dimensions
of the matrix (i.e. the number of bin in the column_dimension and the number
of bins in the row_dimension).

- column_dimension
    - an integer specify how many column bins there are 
- row_dimension
    - an integer specify how many column bins there are 
- column_index
    - column index array
- row_index
    - row index array
- count
    - interaction value measure array
  
run example using:  
> simba3d -r minimial.json

# Why don't you just do a tab or comma delimited format?

Soon! A delimited format is easier to use and more familiar, but it is not as flexible 
as json. Also, with json, you could additionally store other information 
(such as a genomic bin index, a descriptive name for the data, e.c.t.) without 
breaking anything. 


# Is JSON some python thing?

Actually, it is a javascript thing. It happens to be a very minimalistic data
format that is human readable and flexible. Python has libraries for reading
json, and so does every other programing language.

Similiar data formats include XML, which is more verbose, but can do a little 
more than JSON.
