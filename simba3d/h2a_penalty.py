"""



@author: Michael Rosenthal
"""

# import the gradients


import numpy as np

from simba3d.pairwise_computations import triu_ind_to_ij,triu_ij_to_ind


def gradient_h2a(
				penalty_weight,
				length,
				pairwise_difference_x,
				pairwise_difference_y,
				pairwise_difference_z,
				pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient,
				series
				):
	'''
	pairwise repulsion penalty
	'''
	penalty_weight/=(length-2)
	series[0]=0
	for ii in range(1,length-1):
		adj_ind0=int(((ii+1)*(ii))/2.0 -1)# get the adjacent indexes
		adj_ind1=int(((ii+2)*(ii+1))/2.0 -1)# get the adjacent indexes
		dot =(pairwise_difference_x[adj_ind0]*pairwise_difference_x[adj_ind1])
		dot +=(pairwise_difference_y[adj_ind0]*pairwise_difference_y[adj_ind1])
		dot +=(pairwise_difference_z[adj_ind0]*pairwise_difference_z[adj_ind1])
		series[0]-=dot/(pairwise_distance[adj_ind0]*pairwise_distance[adj_ind1])

		if ii==1:
			# compute the first gradient term
			coef0 =-pairwise_difference_x[adj_ind0]*pairwise_difference_x[adj_ind1]
			coef0+=-pairwise_difference_y[adj_ind0]*pairwise_difference_y[adj_ind1]
			coef0+=-pairwise_difference_z[adj_ind0]*pairwise_difference_z[adj_ind1]
			coef0/=pow(pairwise_distance[adj_ind0],3.0)*pairwise_distance[adj_ind1]		

			x_gradient[0]=-penalty_weight*pairwise_difference_x[adj_ind0]*coef0
			y_gradient[0]=-penalty_weight*pairwise_difference_y[adj_ind0]*coef0
			z_gradient[0]=-penalty_weight*pairwise_difference_z[adj_ind0]*coef0
			coef=(pairwise_distance[adj_ind0]*pairwise_distance[adj_ind1])

			x_gradient[0]-=penalty_weight*pairwise_difference_x[adj_ind1]/coef
			y_gradient[0]-=penalty_weight*pairwise_difference_y[adj_ind1]/coef
			z_gradient[0]-=penalty_weight*pairwise_difference_z[adj_ind1]/coef
		elif ii==2:		
			# compute the second gradient term
			# np.sum(pairwise_difference_x[adj_ind0]*pairwise_difference_x[adj_ind1])/(pow(w[0],3.0)*w[1]))
			adj_ind0=int(((2)*(1))/2.0 -1)# get the adjacent indexes
			adj_ind1=int(((3)*(2))/2.0 -1)# get the adjacent indexes
			adj_ind2=int(((4)*(3))/2.0 -1)# get the adjacent indexes
			coef0 =pairwise_difference_x[adj_ind0]*pairwise_difference_x[adj_ind1]
			coef0+=pairwise_difference_y[adj_ind0]*pairwise_difference_y[adj_ind1]
			coef0+=pairwise_difference_z[adj_ind0]*pairwise_difference_z[adj_ind1]
			coef1 =(pow(pairwise_distance[adj_ind0],3.0)*pairwise_distance[adj_ind1])	
			coef2 =1/(pairwise_distance[adj_ind0]*pairwise_distance[adj_ind1])
			coef  =(coef0/coef1+coef2)
			x_gradient[1]-=penalty_weight*pairwise_difference_x[adj_ind0]*coef
			y_gradient[1]-=penalty_weight*pairwise_difference_y[adj_ind0]*coef
			z_gradient[1]-=penalty_weight*pairwise_difference_z[adj_ind0]*coef
	
			coef3 =(pairwise_distance[adj_ind0]*pow(pairwise_distance[adj_ind1],3))	
			coef4 =pairwise_difference_x[adj_ind1]*pairwise_difference_x[adj_ind2]
			coef4+=pairwise_difference_y[adj_ind1]*pairwise_difference_y[adj_ind2]
			coef4+=pairwise_difference_z[adj_ind1]*pairwise_difference_z[adj_ind2]
			coef5 =(pow(pairwise_distance[adj_ind1],3.0)*pairwise_distance[adj_ind2])
			coef  =(-coef2-coef0/coef3- coef4/coef5)
			x_gradient[1]+=-penalty_weight*pairwise_difference_x[adj_ind1]*coef
			y_gradient[1]+=-penalty_weight*pairwise_difference_y[adj_ind1]*coef
			z_gradient[1]+=-penalty_weight*pairwise_difference_z[adj_ind1]*coef

			coef6  =1/(pairwise_distance[adj_ind1]*pairwise_distance[adj_ind2])
			x_gradient[1]+=-penalty_weight*pairwise_difference_x[adj_ind2]*coef6
			y_gradient[1]+=-penalty_weight*pairwise_difference_y[adj_ind2]*coef6
			z_gradient[1]+=-penalty_weight*pairwise_difference_z[adj_ind2]*coef6
		else:
			adj_indn3=int(((ii+-1)*(ii-2))/2.0 -1)# get the adjacent indexes
			adj_indn2=int(((ii+0)*(ii-1))/2.0 -1)# get the adjacent indexes
			adj_indn1=int(((ii+1)*(ii+0))/2.0 -1)# get the adjacent indexes
			adj_ind0=int(((ii+2)*(ii+1))/2.0 -1)# get the adjacent indexes
			coef0  =1/(pairwise_distance[adj_indn3]*pairwise_distance[adj_indn2])
			x_gradient[ii-1]+=penalty_weight*pairwise_difference_x[adj_indn3]*coef0
			y_gradient[ii-1]+=penalty_weight*pairwise_difference_y[adj_indn3]*coef0
			z_gradient[ii-1]+=penalty_weight*pairwise_difference_z[adj_indn3]*coef0
			coef1 =pairwise_difference_x[adj_indn3]*pairwise_difference_x[adj_indn2]
			coef1+=pairwise_difference_y[adj_indn3]*pairwise_difference_y[adj_indn2]
			coef1+=pairwise_difference_z[adj_indn3]*pairwise_difference_z[adj_indn2]
			coef2 =(pairwise_distance[adj_indn3]*pow(pairwise_distance[adj_indn2],3.0))
			coef3 =1/(pairwise_distance[adj_indn2]*pairwise_distance[adj_indn1])
			coef4 =pairwise_difference_x[adj_indn2]*pairwise_difference_x[adj_indn1]
			coef4+=pairwise_difference_y[adj_indn2]*pairwise_difference_y[adj_indn1]
			coef4+=pairwise_difference_z[adj_indn2]*pairwise_difference_z[adj_indn1]
			coef5=pow(pairwise_distance[adj_indn2],3.0)*pairwise_distance[adj_indn1]
			coef=-(coef1/coef2+coef3+coef4/coef5)
			x_gradient[ii-1]+=penalty_weight*pairwise_difference_x[adj_indn2]*coef
			y_gradient[ii-1]+=penalty_weight*pairwise_difference_y[adj_indn2]*coef
			z_gradient[ii-1]+=penalty_weight*pairwise_difference_z[adj_indn2]*coef
				
			coef6=(pairwise_distance[adj_indn2]*pow(pairwise_distance[adj_indn1],3.0))		
			coef7 =pairwise_difference_x[adj_indn1]*pairwise_difference_x[adj_ind0]
			coef7+=pairwise_difference_y[adj_indn1]*pairwise_difference_y[adj_ind0]
			coef7+=pairwise_difference_z[adj_indn1]*pairwise_difference_z[adj_ind0]		
			coef8=pow(pairwise_distance[adj_indn1],3.0)*pairwise_distance[adj_ind0]
			coef=-(-coef4/coef6-coef7/coef8)
			x_gradient[ii-1]+=penalty_weight*pairwise_difference_x[adj_indn1]*coef
			y_gradient[ii-1]+=penalty_weight*pairwise_difference_y[adj_indn1]*coef
			z_gradient[ii-1]+=penalty_weight*pairwise_difference_z[adj_indn1]*coef	

			coef=-1/(pairwise_distance[adj_indn1]*pairwise_distance[adj_ind0])
			x_gradient[ii-1]+=penalty_weight*pairwise_difference_x[adj_ind0]*coef
			y_gradient[ii-1]+=penalty_weight*pairwise_difference_y[adj_ind0]*coef
			z_gradient[ii-1]+=penalty_weight*pairwise_difference_z[adj_ind0]*coef	
	# compute the next to last gradient 	term
	adj_indn4=int(((length+-2)*(length-3))/2.0 -1)# get the adjacent indexes
	adj_indn3=int(((length+-1)*(length-2))/2.0 -1)# get the adjacent indexes
	adj_indn2=int(((length+-0)*(length-1))/2.0 -1)# get the adjacent indexes	
	coef=1/(pairwise_distance[adj_indn4]*pairwise_distance[adj_indn3])
	x_gradient[length-2]+=penalty_weight*pairwise_difference_x[adj_indn4]*coef
	y_gradient[length-2]+=penalty_weight*pairwise_difference_y[adj_indn4]*coef
	z_gradient[length-2]+=penalty_weight*pairwise_difference_z[adj_indn4]*coef		
	coef1 =pairwise_difference_x[adj_indn4]*pairwise_difference_x[adj_indn3]
	coef1+=pairwise_difference_y[adj_indn4]*pairwise_difference_y[adj_indn3]
	coef1+=pairwise_difference_z[adj_indn4]*pairwise_difference_z[adj_indn3]
	coef2=pow(pairwise_distance[adj_indn3],3.0)*pairwise_distance[adj_indn4]
	coef3=1/(pairwise_distance[adj_indn3]*pairwise_distance[adj_indn2])
	coef4 =pairwise_difference_x[adj_indn3]*pairwise_difference_x[adj_indn2]
	coef4+=pairwise_difference_y[adj_indn3]*pairwise_difference_y[adj_indn2]
	coef4+=pairwise_difference_z[adj_indn3]*pairwise_difference_z[adj_indn2]	
	coef5=(pow(pairwise_distance[adj_indn3],3.0)*pairwise_distance[adj_indn2])
	coef=-(coef1/coef2+coef3+coef4/coef5)
	x_gradient[length-2]+=penalty_weight*pairwise_difference_x[adj_indn3]*coef
	y_gradient[length-2]+=penalty_weight*pairwise_difference_y[adj_indn3]*coef
	z_gradient[length-2]+=penalty_weight*pairwise_difference_z[adj_indn3]*coef
	
	#
	coef6=(pairwise_distance[adj_indn3]*pow(pairwise_distance[adj_indn2],3.0))
	coef=(-coef3- coef4/coef6)
	#coef=(-coef3- coef4/(w[n-3]*pow(w[n-2],3.0)))
	x_gradient[length-2]-=penalty_weight*pairwise_difference_x[adj_indn2]*coef
	y_gradient[length-2]-=penalty_weight*pairwise_difference_y[adj_indn2]*coef
	z_gradient[length-2]-=penalty_weight*pairwise_difference_z[adj_indn2]*coef				
	x_gradient[length-1]+=penalty_weight*pairwise_difference_x[adj_indn3]*coef3
	y_gradient[length-1]+=penalty_weight*pairwise_difference_y[adj_indn3]*coef3
	z_gradient[length-1]+=penalty_weight*pairwise_difference_z[adj_indn3]*coef3	
	coef=coef4/coef6
	x_gradient[length-1]-=penalty_weight*pairwise_difference_x[adj_indn2]*coef
	y_gradient[length-1]-=penalty_weight*pairwise_difference_y[adj_indn2]*coef
	z_gradient[length-1]-=penalty_weight*pairwise_difference_z[adj_indn2]*coef			
	series[0]/=length-2