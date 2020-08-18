from simba3d.gradient_manager import gradient_manager
import numpy as np


# example data
contact_row_index=[0]
contact_col_index=[1]
contact_value=[15]
# set parameters
term_parameters={
				"poisson_weight":1.0,
				"h1_weight":1.0,
				"h2b_weight":1.0
				}

gm=gradient_manager(term_parameters)
gm.set_contact_data(contact_row_index,contact_col_index,contact_value)


x=[0.0, 2.0, 3.0, 0.1,1,2]
y=[0.0, 1.0, 2.0, 0.1,2,3]
z=[0.0, 1.0, 2.0, 0.1,3,4]

gm.compute_gradient(x,y,z)

print(gm.x_gradient)
print(gm.y_gradient)
print(gm.z_gradient)