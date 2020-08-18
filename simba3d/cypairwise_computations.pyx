
import numpy as np
cimport numpy as np

cdef extern from "../simba3d_cpp/src/simba3d.hpp" namespace "SIMBA3D":
	void triu_ij_to_ind(int i , int j,int & ind) except +
	void triu_ind_to_ij(int &i , int & j,int ind) except +
	void run_pairwise_computations( 
                int number_of_pairs,
                double *x,
                double *y,
                double *z,
                int number_of_pairs,
                double *pairwise_difference_x,
                double *pairwise_difference_y,
                double *pairwise_difference_z,
                double *pairwise_distance,
                double *series
                ) except +

def run_pairwise_computations(	int number_of_pairs,
								np.ndarray[double,ndim=1,mode="c"] x not None,
								np.ndarray[double,ndim=1,mode="c"] y not None,
								np.ndarray[double,ndim=1,mode="c"] z not None,
								np.ndarray[double,ndim=1,mode="c"] pairwise_difference_x not None,
								np.ndarray[double,ndim=1,mode="c"] pairwise_difference_y not None,
								np.ndarray[double,ndim=1,mode="c"] pairwise_difference_z not None,		
								np.ndarray[double,ndim=1,mode="c"] series not None
								):
	'''
	Stores shared computations derived from each possible combination of pairs
	'''
	run_pairwise_computations(number_of_pairs,&x[0],&y[0],&z[0], &pairwise_difference_x[0],double &pairwise_difference_y[0],&pairwise_difference_z[0],pairwise_distance[0],&series[0])

def run_adjacent_computations(length,pairwise_distance,series):
	'''
	'''
	series[0]=0.0
	series[1]=0.0
	cdef int ii,jj,ind
	for ii in range(length-1):
		adj_ind=int(((ii+2)*(ii+1))/2.0 -1)# get the adjacent indexes
		series[0]+=pairwise_distance[adj_ind]
		series[1]+=pairwise_distance[adj_ind]*pairwise_distance[adj_ind]
	series[2]=series[0]/(length-1)