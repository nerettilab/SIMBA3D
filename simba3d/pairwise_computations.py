import numpy as np
import scipy.spatial.distance as ssd
from scipy.spatial.distance import pdist
from simba3d.matrixlabish import keyboard
def triu_ij_to_ind(i,j):
  i=np.float(i)+1.0
  j=np.float(j)+1.0
  return int((j-1.0)*(j-2.0)/2.0+i-1.0);


def triu_ind_to_ij(ind):
  ind=np.float(ind)+1.0
  j=int(np.ceil((1.0+np.sqrt(1.0+8.0*(ind)))/2.0))
  i=int(( ind- (j-1.0)*(j-2.0)/2.0))
  i-=1
  j-=1
  return i,j

def diff(x1,x2):
    return x1-x2
def run_pairwise_computations(	number_of_pairs,
								x,
								y,
								z,
								pairwise_difference_x,
								pairwise_difference_y,
								pairwise_difference_z,
								pairwise_distance,
								series
								):
	'''
	Stores shared computations derived from each possible combination of pairs
	'''


	series[0]=0.0
	for ind in range(number_of_pairs):
		ii,jj=triu_ind_to_ij(ind);
		# compute the pairwise differences
		pairwise_difference_x[ind]=x[ii]-x[jj];
		pairwise_difference_y[ind]=y[ii]-y[jj];
		pairwise_difference_z[ind]=z[ii]-z[jj];
		# compute the pairwise pairwise_distance
		pairwise_distance[ind] =pairwise_difference_x[ind]*pairwise_difference_x[ind];
		pairwise_distance[ind]+=pairwise_difference_y[ind]*pairwise_difference_y[ind];
		pairwise_distance[ind]+=pairwise_difference_z[ind]*pairwise_difference_z[ind];
		pairwise_distance[ind]=np.sqrt(pairwise_distance[ind]);
		# sum up the pairwise distance while you are at it
		series[0]+=pairwise_distance[ind];
	series[1]=series[0]/number_of_pairs;

	
def run_adjacent_computations(length,pairwise_distance,series):
	'''
	'''
	series[0]=0.0
	series[1]=0.0
	for ii in range(length-1):
		adj_ind=int(((ii+2)*(ii+1))/2.0 -1)# get the adjacent indexes
		series[0]+=pairwise_distance[adj_ind]
		series[1]+=pairwise_distance[adj_ind]*pairwise_distance[adj_ind]
	series[2]=series[0]/(length-1)

