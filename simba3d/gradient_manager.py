# -*- coding: utf-8 -*-
"""



@author: Michael Rosenthal
"""

# import the gradients
'''
Some gradients have cython code associated with them to speed up the 
computation. Try and import the cython and import a python if that fails. 
'''

import numpy as np
#import time
#import scipy.spatial.distance as ssd
#from scipy.spatial.distance import pdist
#from scipy.sparse import csr_matrix
#from scipy.optimize import minimize,check_grad
#from scipy.misc import comb
import scipy.sparse as sp

from simba3d.pairwise_computations import run_pairwise_computations,run_adjacent_computations
from simba3d.gradient_poisson import gradient_poisson
from simba3d.h1_penalty import gradient_h1    
from simba3d.h2a_penalty import gradient_h2a   
from simba3d.h2b_penalty import run_h2b_computations,gradient_h2b   
from simba3d.h2c_penalty import run_h2c_computations,gradient_h2c  

class gradient_manager():
	'''
	Manages the gradients 

	As simba3d has developed, we have gone through many penalty functions. To 
	ease the burden of maintaining these functions, a  class object 
	has been made to efficiciently manage gradient computations.

	The reason this is needed is because the gradient computations need to be
	very efficient. They are called thousands of times within an optimization. 
	This requires efficient coding and resource utlization. Several functions
	use the same types of computations. Rather than computing these values
	multiple times, it is better to store them within the object class and then
	utilize them as required for a particular setting.

	There is also the use of c and cython code. As a policy, we first create
	a python version of the gradient, and then create a niave cython implementation.
	Afterwards, the cython code is further refined to boost performance. This requires
	copilation of source code. If the cython code fails to import, there is a python
	implementation which can be loaded natively with a minimal number of depencies.

	'''
	def __init__(self,term_parameters=None):
		# initialize the stored inputs, calculations, and parameters
		# set default values for parameters here (set to None to error out when misspecified)
		self.contact_row_index=None
		self.contact_col_index=None
		self.contact_value=None
		self.contact_sum=None
		self.lenght=None
		self.number_of_pairs=None
		self.pairwise_difference_x=list()
		self.pairwise_difference_y=list()
		self.pairwise_difference_z=list()
		self.pairwise_distance=list()
		self.sum_pairwise_distance=None
		self.average_pairwise_distance=None
		self.sum_adjacent_distance=None
		self.sum_of_squared_adjacent_distance=None
		self.average_adjacent_distance=None
		self.poisson_a_parameter=-3.0
		self.poisson_b_parameter=1.0
		self.poisson_loglikelihood=None
		self.h1_uniform_spacing_value=None
		self.h2a_value=None
		self.h2b_alpha=10
		self.h2b_radius=1
		self.h2b_unwieghted_value=None
		self.h2b_F=None
		self.h2b_H=None
		self.h2b_G=None
		self.h2c_radius=1.0
		self.h2c_unwieghted_value=None
		self.x_gradient=None
		self.y_gradient=None
		self.z_gradient=None
		self.e=None		
		# a list of functions is created upon initialization to eliminate
		# the need for if statements or switches during the optimization
		self.functions =[]
		# For each iteration, there are calculations that are used in multiple gradients. 
		# This list tells the gradient manager which computations it needs to update prior
		# to computing the gradients.
		self.shared_computations=[]
		pairwise_computations_needed=False
		adjacent_computations_needed=False
		h2b_computations_needed=False
		h2c_computations_needed=False
		if term_parameters is None:
			term_parameters={}
		if "poisson_weight" in term_parameters:
			if term_parameters["poisson_weight"]>0:
				self.functions.append(["poisson",term_parameters["poisson_weight"],self.gradient_data])
				pairwise_computations_needed=True
			if "poisson_parameter_a" in term_parameters:
				 self.poisson_a_parameter=term_parameters["poisson_parameter_a"]
			if "poisson_parameter_b" in term_parameters:
				 self.poisson_b_parameter=term_parameters["poisson_parameter_b"]
		if "h1_weight" in term_parameters:
			if term_parameters["h1_weight"]>0:
				self.functions.append(["h1",term_parameters["h1_weight"],self.gradient_h1])
				pairwise_computations_needed=True
				adjacent_computations_needed=True
		if "h2a_weight" in term_parameters:
			if term_parameters["h2a_weight"]>0:
				self.functions.append(["h2a",term_parameters["h2a_weight"],self.gradient_h2a])
				pairwise_computations_needed=True
				adjacent_computations_needed=True
			if "h2a_parameter_radius" in term_parameters:
				self.h2c_radius=term_parameters["h2a_parameter_radius"]
		if "h2b_weight" in term_parameters:
			if term_parameters["h2b_weight"]>0:
				self.functions.append(["h2b",term_parameters["h2b_weight"],self.gradient_h2b])
				pairwise_computations_needed=True
				adjacent_computations_needed=True
				h2b_computations_needed=True
		if "h2c_weight" in term_parameters:
			if term_parameters["h2c_weight"]>0:
				self.functions.append(["h2c",term_parameters["h2c_weight"],self.gradient_h2c])
				pairwise_computations_needed=True
				adjacent_computations_needed=True
				h2c_computations_needed=True
		# append the shared computation tasks to the list of shared computations
		if pairwise_computations_needed:
			self.shared_computations.append(self.run_pairwise_computations)
		if adjacent_computations_needed:
			self.shared_computations.append(self.run_adjacent_computations)
		if h2b_computations_needed:
			self.shared_computations.append(self.run_h2b_computations)
		if h2c_computations_needed:
			self.shared_computations.append(self.run_h2c_computations)

	def set_contact_data(self,contact_row_index,contact_col_index,contact_value):
		'''
		'''
		self.contact_row_index=contact_row_index
		self.contact_col_index=contact_col_index
		self.contact_value=contact_value
		self.contact_sum=sum(contact_value)
	def add_population_contact_data(self,weight,contact_row_index,contact_col_index,contact_value):
		'''
		'''
		C=sp.coo_matrix((self.ontact_value,(self.contact_row_index,self.contact_col_index)),shape=(self.length,self.length))
		Cpop=sp.coo_matrix((contact_value,(contact_row_index,contact_col_index)),shape=(self.length,self.length))
		C=C+weight*Cpop
		self.contact_row_index=C.row.tolist()
		self.contact_col_index=C.col.tolist()
		self.contact_value=C.data.tolist()
		self.contact_sum=sum(self.contact_value)
	def run_pairwise_computations(self):
		'''
		Stores shared computations derived from each possible combination of pairs
		'''
		# pass in mutable objects
		self.pairwise_difference_x=[0]*self.number_of_pairs
		self.pairwise_difference_y=[0]*self.number_of_pairs
		self.pairwise_difference_z=[0]*self.number_of_pairs
		self.pairwise_distance=[0]*self.number_of_pairs
		series=[0]*2 
		# There are several scalar values computed. Rather than outputing them in the return
		# I opted to pass in mutable objects so that the memory is already allocated
		# for them and I can assume the memory is allocated when I make the c code.
		run_pairwise_computations(	
								self.number_of_pairs,
								self.x,
								self.y,
								self.z,
								self.pairwise_difference_x,
								self.pairwise_difference_y,
								self.pairwise_difference_z,
								self.pairwise_distance,
								series
								)
		self.sum_pairwise_distance = series[0]
		self.average_pairwise_distance= series[1]

	def run_adjacent_computations(self):
		'''
		Stores shared computations derived from each adjacent pair of points

		Depends on computations made during run_pairwise_computations
		'''
		series=[0]*4
		run_adjacent_computations(self.length,self.pairwise_distance,series)
		self.sum_adjacent_distance = series[0]
		self.sum_of_squared_adjacent_distance= series[1]
		self.average_adjacent_distance= series[2]
	def run_h2b_computations(self):
		'''

		Depends on computations made during run_adjacent_computations and run_pairwise_computations
		'''
		series=[0]*1
		self.h2b_F=[0]*	self.number_of_pairs
		self.h2b_G=[0]*	self.number_of_pairs
		self.h2b_H=[0]*	self.number_of_pairs
		run_h2b_computations(
							self.length,
							self.number_of_pairs,
							self.h2b_alpha,
							self.average_adjacent_distance,
							self.pairwise_distance,
							self.h2b_F,
							self.h2b_G,
							self.h2b_H,
							series
							)
		self.h2b_unwieghted_value=series[0]
	def run_h2c_computations(self):
		'''

		Depends on computations made during run_adjacent_computations and run_pairwise_computations
		'''
		series=[0]*1
		self.h2b_F=[0]*	self.number_of_pairs
		run_h2c_computations(
							self.length,
							self.number_of_pairs,
							self.h2c_radius,
							self.average_adjacent_distance,
							self.pairwise_distance,
							self.h2b_F,
							series
							)
		self.h2c_unwieghted_value=series[0]
	def gradient_data(self,penalty_weight,x_gradient,y_gradient,z_gradient):
		'''
		'''
		series=[0]*1
		gradient_poisson(
						penalty_weight,
						self.contact_row_index,
						self.contact_col_index,
						self.contact_value,
						self.poisson_a_parameter,
						self.poisson_b_parameter,
						self.contact_sum,
						self.pairwise_difference_x,
						self.pairwise_difference_y,
						self.pairwise_difference_z,
						self.pairwise_distance,
						x_gradient,
						y_gradient,
						z_gradient,
						series
						)
		self.poisson_loglikelihood  =penalty_weight*series[0]
		return  self.poisson_loglikelihood

	def gradient_h1(self,penalty_weight,x_gradient,y_gradient,z_gradient):	
		'''
		A uniform spacing penalty which will force the nodes to be spaces uniformly.

		Depends on computations made during run_adjacent_computations and run_pairwise_computations
		'''
		gradient_h1(				
				penalty_weight,
				self.length,
				self.sum_adjacent_distance,
				self.sum_of_squared_adjacent_distance,
				self.pairwise_difference_x,
				self.pairwise_difference_y,
				self.pairwise_difference_z,
				self.pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient
					)
		self.h1_uniform_spacing_value = penalty_weight*((self.length-1)*self.sum_of_squared_adjacent_distance/(self.sum_adjacent_distance*self.sum_adjacent_distance) -1)
		return self.h1_uniform_spacing_value
	def gradient_h2a(self,penalty_weight,x_gradient,y_gradient,z_gradient):	
		'''
		'''
		series=[0]*1
		gradient_h2a(
				penalty_weight,
				self.length,
				self.pairwise_difference_x,
				self.pairwise_difference_y,
				self.pairwise_difference_z,
				self.pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient,
				series
				)
		self.h2a_value=penalty_weight*series[0];
		return self.h2a_value		
	def gradient_h2b(self,penalty_weight,x_gradient,y_gradient,z_gradient):	
		'''
		'''
		gradient_h2b(
				penalty_weight,
				self.length,
				self.number_of_pairs,
				self.h2b_alpha,
				self.h2b_F,
				self.h2b_G,
				self.h2b_H,
				self.average_adjacent_distance,
				self.pairwise_difference_x,
				self.pairwise_difference_y,
				self.pairwise_difference_z,
				self.pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient
				)
		return penalty_weight*self.h2b_unwieghted_value
	def gradient_h2c(self,penalty_weight,x_gradient,y_gradient,z_gradient):	
		'''
		'''
		gradient_h2c(
				penalty_weight,
				self.length,
				self.number_of_pairs,
				self.h2b_alpha,
				self.h2b_F,
				self.h2c_radius,
				self.average_adjacent_distance,
				self.pairwise_difference_x,
				self.pairwise_difference_y,
				self.pairwise_difference_z,
				self.pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient
				)
		return penalty_weight*self.h2c_unwieghted_value
	def compute_gradient(self,x,y,z,x_gradient,y_gradient,z_gradient):
		'''
		Compute the energy and gradient of the functions for a given input curve
		'''
		self.length=len(x)
		self.number_of_pairs=int(self.length*(self.length-1)/2);
		self.x=list(x)
		self.y=list(y)
		self.z=list(z)
		# Loop through the list of shared calculations
		for shared_computation in self.shared_computations:
			shared_computation()

		# initialize the gradient
		self.x_gradient=[0]*self.length
		self.y_gradient=[0]*self.length
		self.z_gradient=[0]*self.length
		# initialize total energy
		self.e=0.0
		# loop through the non-trivial terms
		for functional in self.functions:
			self.e+=functional[2](functional[1],x_gradient,y_gradient,z_gradient)