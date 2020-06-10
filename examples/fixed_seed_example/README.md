You can set the seed, but I think it is better to set the initialized curve if
you want your results to be reproducable.

To go through this example, run the following command twice.
> simba3d -i tasks/*

Then to compare them run the following command:
> simba3d-result-disp -p fixed -i *.json

You should get two identical results with different names and computation times.
