try:
	from simba3d.cysimba3d import gradient_manager
except:
	print('Failed to import cython class (you probably do not have it compiled for your system).')
	print('Loading python version (which is slower and does not use multi-threading)')
	from simba3d.gradient_manager import gradient_manager

import numpy as np

if __name__ == "__main__":
	# example data
	# example data
	contact_row_index=[0]
	contact_col_index=[1]
	contact_value=[15]
	# set parameters
	term_parameters={
					"poisson_weight":1.0,
					"h1_weight":1.0,
					"h2c_weight":1.0
					}

	gm=gradient_manager(term_parameters)
	gm.set_contact_data(contact_row_index,contact_col_index,contact_value)


	x=[0.0, 2.0, 3.0, 0.1,1,2]
	y=[0.0, 1.0, 2.0, 0.1,2,3]
	z=[0.0, 1.0, 2.0, 0.1,3,4]

	# allocate memory for the numpy array
	x_gradient=np.zeros(len(x))
	y_gradient=np.zeros(len(x))
	z_gradient=np.zeros(len(x))
	gm.compute_gradient(x,y,z,x_gradient,y_gradient,z_gradient)
	print(x_gradient)
	print(y_gradient)
	print(z_gradient)