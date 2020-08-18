"""



@author: Michael Rosenthal
"""

# import the gradients


import numpy as np

from simba3d.pairwise_computations import triu_ind_to_ij,triu_ij_to_ind
def run_h2b_computations(
							length,
							number_of_pairs,
							h2b_alpha,
							average_adjacent_distance,
							pairwise_distance,
							h2b_F,
							h2b_G,
							h2b_H,
							series
							):
	'''
	'''
	series[0]=0.0
	for ind in range(number_of_pairs):
		ii,jj=triu_ind_to_ij(ind)
		u_jk= pairwise_distance[ind]/average_adjacent_distance
		h2b_F[ind]=1.0-u_jk
		h2b_G[ind]=np.exp(h2b_alpha*(u_jk-1.0));
		h2b_H[ind]=1.0/(1.0+h2b_G[ind]);
		series[0]+=h2b_F[ind]*h2b_F[ind]*h2b_F[ind]*h2b_F[ind]*h2b_H[ind];
	series[0]*=2.0/np.float( (length*(length-1.0)));

def gradient_h2b(
				penalty_weight,
				length,
				number_of_pairs,
				h2b_alpha,
				h2b_F,
				h2b_G,
				h2b_H,
				average_adjacent_distance,
				pairwise_difference_x,
				pairwise_difference_y,
				pairwise_difference_z,
				pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient
				):
	'''
	pairwise repulsion penalty
	'''
	for ii in range(length):
		# partial derivative of mu
		if ii==0:
			ind2=triu_ij_to_ind(0,1);
			partial_mu_xi_x=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
		elif (ii==length-1):
			ind2=triu_ij_to_ind(ii-1,ii);
			partial_mu_xi_x=-pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=-pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=-pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
		else:
			ind2=triu_ij_to_ind(ii,ii+1);
			partial_mu_xi_x=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
			ind2=triu_ij_to_ind(ii-1,ii);
			partial_mu_xi_x-=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y-=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z-=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);		
		for ind in range(number_of_pairs):
			jj,kk=triu_ind_to_ij(ind);
			# partial derivative of pairwise_distance
			partial_djk_xi_x=0.0;
			partial_djk_xi_y=0.0;
			partial_djk_xi_z=0.0;
			if (jj==ii):
				partial_djk_xi_x=pairwise_difference_x[ind]/pairwise_distance[ind];
				partial_djk_xi_y=pairwise_difference_y[ind]/pairwise_distance[ind];
				partial_djk_xi_z=pairwise_difference_z[ind]/pairwise_distance[ind];
			elif (kk==ii):
				partial_djk_xi_x=-pairwise_difference_x[ind]/pairwise_distance[ind];
				partial_djk_xi_y=-pairwise_difference_y[ind]/pairwise_distance[ind];
				partial_djk_xi_z=-pairwise_difference_z[ind]/pairwise_distance[ind];
			# partial derivative of ujk
			partial_ujk_xi_x=(partial_djk_xi_x+(h2b_F[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
			partial_ujk_xi_y=(partial_djk_xi_y+(h2b_F[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
			partial_ujk_xi_z=(partial_djk_xi_z+(h2b_F[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
			# compute gradient
			scale=-penalty_weight*2/(length*(length-1));
			term=scale*h2b_F[ind]*h2b_F[ind]*h2b_F[ind]*h2b_H[ind]*(4.0+h2b_alpha*h2b_F[ind]*h2b_G[ind]*h2b_H[ind]);
			x_gradient[ii]+=term*partial_ujk_xi_x;
			y_gradient[ii]+=term*partial_ujk_xi_y;
			z_gradient[ii]+=term*partial_ujk_xi_z;