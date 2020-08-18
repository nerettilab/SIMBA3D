/**
* @file simba3d.cpp
* @brief a brief description for  the template name cpp file 
*
* If you want, you can say more about the contents of the header. 
*
*
*/

#include "simba3d.hpp"

#include "pthread.h"
#include  <cmath> 
//#include  <stdio.h> 
//#include "mex.h"
#include <ctime>

void *split_compute_pairwise_computations_helper(
  void *args
  ){
  	thread_struct *info = (thread_struct *) args;
    info->instance->split_compute_pairwise_computations(info->subtask_id);
	return args;
}

void *split_compute_poisson_loglikelihood_gradient_helper(
  void *args
  ){
  	thread_struct *info = (thread_struct *) args;
    info->instance->split_compute_poisson_loglikelihood_gradient(info->subtask_id,info->grad_x,info->grad_y,info->grad_z);
	return args;
}
void *split_compute_h2c_gradient_helper(
  void *args
  ){
  	thread_struct *info = (thread_struct *) args;
    info->instance->split_compute_h2c_gradient(info->subtask_id,info->grad_x,info->grad_y,info->grad_z);
	return args;
}
void SIMBA3D::triu_ij_to_ind(int i , int j,int & ind){
  i++;
  j++;
  ind=(int) ( ((double) j-1.0)*((double) j-2.0)/2.0+i-1);
}

void SIMBA3D::triu_ind_to_ij(int &i , int & j,int ind){
  ind++;
  j=std::ceil((1.0+std::sqrt(1.0+8.0*((double) ind)))/2.0);
  i= (int) ((double) ind- (j-1.0)*(j-2.0)/2.0);
  i--;
  j--;
}

void  SIMBA3D::compute_pairwise_difference(
	const std::vector<double> & input_x,
	/**< \brief x coordinate*/
	const std::vector<double> & input_y,
	/**< \brief y coordinate*/
	const std::vector<double> & input_z,
	/**< \brief z coordinate*/
	std::vector<double> & x_differences,
	/**< \brief upper triangular differences*/
	std::vector<double> & y_differences,
	/**< \brief upper triangular differences*/
	std::vector<double> & z_differences,
	/**< \brief upper triangular differences*/
	std::vector<double> & distances
	/**< \brief upper triangular differences*/
){
	int n=(int) input_x.size();
	n=std::min(n,(int) input_y.size());
	n=std::min(n,(int) input_z.size());
	int m=n*(n-1)/2;// dimension of upper triangular matrix

	if (x_differences.size()!=(unsigned) m){
		x_differences.clear();
		x_differences.reserve(m);
		x_differences.resize(m);
	}
	if (y_differences.size()!=(unsigned) m){
		y_differences.clear();
		y_differences.reserve(m);
		y_differences.resize(m);
	}
	if (z_differences.size()!=(unsigned) m){
		z_differences.clear();
		z_differences.reserve(m);
		z_differences.resize(m);
	}
	if (distances.size()!=(unsigned) m){
		distances.clear();
		distances.reserve(m);
		distances.resize(m);
	}
	int ind,ii,jj;
	for (ind=0;ind<m;ind++){
		SIMBA3D::triu_ind_to_ij(ii ,jj,ind);
		x_differences[ind]=input_x[ii]-input_x[jj];
		y_differences[ind]=input_y[ii]-input_y[jj];
		z_differences[ind]=input_z[ii]-input_z[jj];
		distances[ind]=x_differences[ind]*x_differences[ind];

		distances[ind]+=y_differences[ind]*y_differences[ind];
		distances[ind]+=z_differences[ind]*z_differences[ind];
		distances[ind]=sqrt(distances[ind]);
	}
}
//.................................................................................

SIMBA3D::simba3d::simba3d(){
	penalty_inventory_length=5;// the number of penalty terms
	initialize_penalty_weights();
	length=0;
	poisson_a_parameter=-3.0;
    poisson_b_parameter=1.0;
    alpha_h2b=10.0;
    h2c_radius=1.0;
  
}
//.................................................................................

SIMBA3D::simba3d::~simba3d() {
  //  deconstruct the dynamically allocated variables (any variable you used new)
}
//.................................................................................
void SIMBA3D::simba3d::initialize_penalty_weights(){
	penalty_weights.clear();
	
	for (int ii=0; ii <penalty_inventory_length ; ii++){
		penalty_weights.push_back(0.0);
		penalty_weights.push_back(0.0);
		penalty_weights.push_back(0.0);
		penalty_weights.push_back(0.0);
		penalty_weights.push_back(0.0);		
	}
	penalty_weights[0]=1.0;


}
//.................................................................................
void SIMBA3D::simba3d::set_penalty_weights(
      const std::vector<double> & input_penalty_weights
){
	penalty_weights.clear();
	for (unsigned ii=0;ii<input_penalty_weights.size();ii++){
		penalty_weights.push_back(input_penalty_weights[ii]);  
	}
	// if the user did not provide enough penalty weight coefficients,
	// then make the rest zero
	while (penalty_weights.size() < (unsigned) penalty_inventory_length){
		penalty_weights.push_back(0.0);
	}	
}
void SIMBA3D::simba3d::set_penalty_weights(
   int input_number_of_penalties,
   double * input_penalty_weights
   ){
	penalty_weights.clear();
	
	penalty_weights.insert(penalty_weights.begin(),input_penalty_weights,input_penalty_weights+input_number_of_penalties);
	//penalty_weights.resize(input_number_of_penalties);
	// if the user did not provide enough penalty weight coefficients,
	// then make the rest zero
	while (penalty_weights.size() < (unsigned) penalty_inventory_length){
		penalty_weights.push_back(0.0);
	}	
}
void SIMBA3D::simba3d::set_penalty_parameters(
   int input_number_of_penalties,
   double * input_penalty_parameters
   ){
	   if (input_number_of_penalties>0){
		   poisson_a_parameter=input_penalty_parameters[0];
	   }
	   if (input_number_of_penalties>1){
		   poisson_b_parameter=input_penalty_parameters[1];
	   }
	   if (input_number_of_penalties>2){
		   alpha_h2b=input_penalty_parameters[2];
	   }
	   if (input_number_of_penalties>3){
		   h2c_radius=input_penalty_parameters[3];
	   }
}
//.................................................................................
void SIMBA3D::simba3d::set_curve(
    const std::vector<double> & input_x,
    /**< \brief x coordinate*/
    const std::vector<double> & input_y,
    /**< \brief y coordinate*/
    const std::vector<double> & input_z
    /**< \brief z coordinate*/
){
	length=(int) input_x.size();
	length=std::min(length,(int) input_y.size());
	length=std::min(length,(int) input_z.size());
	number_of_pairs=length*(length-1)/2;
	x.clear();
	y.clear();
	z.clear();
	x=input_x;
	y=input_y;
	z=input_z;
}
void SIMBA3D::simba3d::set_curve(
	int input_length,
    double * input_x,
    /**< \brief x coordinate*/
    double * input_y,
    /**< \brief y coordinate*/
    double * input_z
    /**< \brief z coordinate*/
	){
	length=input_length;
	number_of_pairs=length*(length-1)/2;
	x.clear();
	y.clear();
	z.clear();
	x.insert(x.begin(),input_x,input_x+length);
	y.insert(y.begin(),input_y,input_y+length);
	z.insert(z.begin(),input_z,input_z+length);

}
void SIMBA3D::simba3d:: set_pairwise(
	int input_length,
    double * input_x_differences,
    /**< \brief upper triangular differences*/
    double * input_y_differences,
    /**< \brief upper triangular differences*/
    double * input_z_differences,
   /**< \brief upper triangular differences*/   
    double * input_distances
   /**< \brief upper triangular distances*/
  ){
	length=input_length;
	number_of_pairs=length*(length-1)/2;
	x.clear();
	y.clear();
	z.clear();
	x_differences.clear();
	y_differences.clear();
	z_differences.clear();
	distances.clear();
	x_differences.insert(x_differences.begin(),input_x_differences,input_x_differences+length);
	y_differences.insert(y_differences.begin(),input_y_differences,input_y_differences+length);
	z_differences.insert(z_differences.begin(),input_z_differences,input_z_differences+length);
	distances.insert(distances.begin(),input_distances,input_distances+length);
}
void SIMBA3D::simba3d::set_contact(
    const std::vector<double> & row_index,
    /**< \brief x coordinate*/
    const std::vector<double> & column_index,
    /**< \brief y coordinate*/
    const std::vector<double> & value
    /**< \brief z coordinate*/
  ){
	number_of_contact_pairs=(int) row_index.size();
	number_of_contact_pairs=std::min(number_of_contact_pairs,(int) column_index.size());
	number_of_contact_pairs=std::min(number_of_contact_pairs,(int) value.size());
	contact_sum=0.0;
	if (contact_row_index.size()!= (unsigned) number_of_contact_pairs){
		// if the dimension does not match, clear out the memory
		contact_row_index.clear();
		contact_column_index.clear();
		contact_value.clear();
		contact_row_index.reserve(number_of_contact_pairs);
		contact_column_index.reserve(number_of_contact_pairs);
		contact_value.reserve(number_of_contact_pairs);
		contact_row_index.resize(number_of_contact_pairs);
		contact_column_index.resize(number_of_contact_pairs);
		contact_value.resize(number_of_contact_pairs);
	}
	for (int ii=0; ii<number_of_contact_pairs;ii++){
		if (row_index[ii]<column_index[ii]){
			contact_row_index[ii]=row_index[ii];
			contact_column_index[ii]=column_index[ii];
			contact_value[ii]=value[ii];
			contact_sum+=value[ii];
		} else if (row_index[ii]>column_index[ii]){
			contact_row_index[ii]=column_index[ii];
			contact_column_index[ii]=row_index[ii];
			contact_value[ii]=value[ii];
			contact_sum+=value[ii];
		} else {
			number_of_contact_pairs--;
		}
	}

}
void SIMBA3D::simba3d::set_contact(
    int input_number_of_contacts,
    double * row_index,
    double * column_index,
    double * value
  ){
	number_of_contact_pairs=input_number_of_contacts;
	contact_row_index.clear();
	contact_column_index.clear();
	contact_value.clear();
	contact_row_index.insert(contact_row_index.begin(),row_index,row_index+number_of_contact_pairs);
	contact_column_index.insert(contact_column_index.begin(),column_index,column_index+number_of_contact_pairs);
	contact_value.insert(contact_value.begin(),value,value+number_of_contact_pairs);
	contact_sum=0.0;
	for (int ii=0; ii<number_of_contact_pairs;ii++){
		if (row_index[ii]<column_index[ii]){
			contact_row_index[ii]=row_index[ii];
			contact_column_index[ii]=column_index[ii];
			contact_value[ii]=value[ii];
			contact_sum+=value[ii];
		} else if (row_index[ii]>column_index[ii]){
			contact_row_index[ii]=column_index[ii];
			contact_column_index[ii]=row_index[ii];
			contact_value[ii]=value[ii];
			contact_sum+=value[ii];
		} else {
			number_of_contact_pairs--;
		}
	}

}
int SIMBA3D::simba3d::get_contact_pairs(unsigned ind,int & r, int &c, int &v){
	if (ind <(unsigned) number_of_contact_pairs) {
		r=contact_row_index[ind];
		c=contact_column_index[ind];
		v=contact_value[ind];
		return 0;
	} else {
		return -1;
	}
}
//.................................................................................
void SIMBA3D::simba3d::compute_pairwise_computations()
{
	//initialize_pairwise_computations();

	int ind,ii,jj;
	sum_pairwise_distance=0.0;


	for (ind=0;ind<number_of_pairs;ind++){
		SIMBA3D::triu_ind_to_ij(ii ,jj,ind);
		x_differences[ind]=x[ii]-x[jj];
		y_differences[ind]=y[ii]-y[jj];
		z_differences[ind]=z[ii]-z[jj];
		distances[ind]=x_differences[ind]*x_differences[ind];
		distances[ind]+=y_differences[ind]*y_differences[ind];
		distances[ind]+=z_differences[ind]*z_differences[ind];
		distances[ind]=sqrt(distances[ind]);
		sum_pairwise_distance+=distances[ind];
	}
	average_pairwise_distance=sum_pairwise_distance/number_of_pairs;

}
void SIMBA3D::simba3d::initialize_pairwise_computations(){
	if (x_differences.size()!= (unsigned) number_of_pairs){
		x_differences.clear();
		x_differences.reserve(number_of_pairs);
		x_differences.resize(number_of_pairs);
	}
	if (y_differences.size()!=(unsigned) number_of_pairs){
		y_differences.clear();
		y_differences.reserve(number_of_pairs);
		y_differences.resize(number_of_pairs);
	}
	if (z_differences.size()!=(unsigned) number_of_pairs){
		z_differences.clear();
		z_differences.reserve(number_of_pairs);
		z_differences.resize(number_of_pairs);
	}
	if (distances.size()!=(unsigned) number_of_pairs){
		distances.clear();
		distances.reserve(number_of_pairs);
		distances.resize(number_of_pairs);
	}
}
/*
void SIMBA3D::simba3d::compute_pairwise_computations(int number_of_threads){
  set_number_of_subtasks(number_of_threads);
  // initialize threads
  initialize_pairwise_computations();
  std::thread *pThreads[number_of_threads];
  for (int ii=0;ii<number_of_threads;ii++){
    pThreads[ii]=new std::thread (&SIMBA3D::simba3d::split_compute_pairwise_computations,this,ii);
  }
  // wait for all the threads to finish
  for (int ii=0;ii<number_of_threads;ii++){
     pThreads[ii]->join();
  }
  // wait for all the threads to finish
  for (int ii=0;ii<number_of_threads;ii++){
     delete pThreads[ii];
     pThreads[ii]=0;
  }
  combine_split_pairwise_sums();
}
*/
void SIMBA3D::simba3d::compute_pairwise_computations_0(int number_of_threads){

  set_number_of_subtasks(number_of_threads);
  // initialize threads
  //initialize_pairwise_computations();


  std::vector<pthread_t> pThreads;
  pThreads.clear();
  pThreads.resize(number_of_threads);
  std::vector<thread_struct> tInfo;
  tInfo.clear();
  tInfo.resize(number_of_threads);
  for (int ii=0;ii<number_of_threads;ii++){
	tInfo[ii].instance=this;
  	tInfo[ii].subtask_id=ii;
    pthread_create(&pThreads[ii],NULL,split_compute_pairwise_computations_helper,(void *) &tInfo[ii]);
  }

  // wait for all the threads to finish
  for (int ii=0;ii<number_of_threads;ii++){
	pthread_join(pThreads[ii],NULL);
  }

  combine_split_pairwise_sums();

}
void SIMBA3D::simba3d::set_number_of_subtasks(int number){

	number_of_subtasks=number;
	split_sum_pairwise_distance.clear();
	split_sum_pairwise_distance.reserve(number_of_subtasks);
	split_sum_pairwise_distance.resize(number_of_subtasks);
	std::fill(split_sum_pairwise_distance.begin(),split_sum_pairwise_distance.end(),0.0);
	//std::cout << number_of_subtasks << " split of the tasks\n";
}

void SIMBA3D::simba3d::split_compute_pairwise_computations(
  int subtask_id
  ){
	//mexPrintf("%d\n",subtask_id);	 
	//std::cout << subtask_id << "\n";
	int ind,ii,jj;
	split_sum_pairwise_distance[subtask_id]=0.0;
	for (ind=subtask_id;ind<number_of_pairs;ind+=number_of_subtasks){
		SIMBA3D::triu_ind_to_ij(ii ,jj,ind);
		x_differences[ind]=x[ii]-x[jj];
		y_differences[ind]=y[ii]-y[jj];
		z_differences[ind]=z[ii]-z[jj];
		distances[ind]=x_differences[ind]*x_differences[ind];
		distances[ind]+=y_differences[ind]*y_differences[ind];
		distances[ind]+=z_differences[ind]*z_differences[ind];
		distances[ind]=sqrt(distances[ind]);;
		split_sum_pairwise_distance[subtask_id]+=distances[ind];
	}
}
void SIMBA3D::simba3d::combine_split_pairwise_sums(){
	sum_pairwise_distance=0.0;
	for (int subtask_id=0;subtask_id<number_of_subtasks;subtask_id++){
		sum_pairwise_distance+=split_sum_pairwise_distance[subtask_id];
	}
	average_pairwise_distance=sum_pairwise_distance/number_of_pairs;
}
//.................................................................................
void SIMBA3D::simba3d::compute_adjacent_computations()
{
	int adj_ind;
	sum_adjacent_distance=0.0;
	sum_of_squared_adjacent_distance=0.0;
	for (int ii=0;ii<length-1;ii++){
		adj_ind=(ii+2)*(ii+1)/2 -1;// get the adjacent indexes
		
		x_differences[adj_ind]=x[ii]-x[ii+1];
		y_differences[adj_ind]=y[ii]-y[ii+1];
		z_differences[adj_ind]=z[ii]-z[ii+1];
		distances[adj_ind]=x_differences[adj_ind]*x_differences[adj_ind];
		distances[adj_ind]+=y_differences[adj_ind]*y_differences[adj_ind];
		distances[adj_ind]+=z_differences[adj_ind]*z_differences[adj_ind];
		distances[adj_ind]=sqrt(distances[adj_ind]);;
		/**/
		sum_adjacent_distance+=distances[adj_ind];
		sum_of_squared_adjacent_distance+=distances[adj_ind]*distances[adj_ind];
		//x_differences[adj_ind];
		//y_differences[adj_ind];
		//z_differences[adj_ind];
		/* if you need to row and column index
		int r,c;
		SIMBA3D::triu_ind_to_ij(r,c,adj_ind);
		*/
	}
	average_adjacent_distance=sum_adjacent_distance/(length-1);
	h1_uniform_spacing_value=(length-1)*sum_of_squared_adjacent_distance/(sum_adjacent_distance*sum_adjacent_distance) -1;
}

//.................................................................................
void SIMBA3D::simba3d::clear_gradient(){
	if (x_gradient.size()!=(unsigned) length){
		x_gradient.clear();
		x_gradient.reserve(length);
		x_gradient.resize(length);
	}
	if (y_gradient.size()!=(unsigned) length){
		y_gradient.clear();
		y_gradient.reserve(length);
		y_gradient.resize(length);
	}
	if (z_gradient.size()!=(unsigned) length){
		z_gradient.clear();
		z_gradient.reserve(length);
		z_gradient.resize(length);
	}
	std::fill(x_gradient.begin(),x_gradient.end(),0.0);
	std::fill(y_gradient.begin(),y_gradient.end(),0.0);
	std::fill(z_gradient.begin(),z_gradient.end(),0.0);
}

void SIMBA3D::simba3d::clear_gradient(double *grad_x,double *grad_y,double *grad_z){
	x_gradient.clear();
	y_gradient.clear();
	z_gradient.clear();
	x_gradient.insert(x_gradient.begin(),grad_x,grad_x+length);
	y_gradient.insert(y_gradient.begin(),grad_y,grad_y+length);
	z_gradient.insert(z_gradient.begin(),grad_z,grad_z+length);


	std::fill(x_gradient.begin(),x_gradient.end(),0.0);
	std::fill(y_gradient.begin(),y_gradient.end(),0.0);
	std::fill(z_gradient.begin(),z_gradient.end(),0.0);

}
void SIMBA3D::simba3d::print_gradient(){
	std::cout <<"gradient"<<"\n";
	for (unsigned ii=0;ii<x_gradient.size();ii++){
		std::cout << x_gradient[ii] << " " << y_gradient[ii] << " " << z_gradient[ii] << "\n";
	}
}

void SIMBA3D::simba3d::compute_poisson_loglikelihood(){
	if (penalty_weights[0]>0){
	poisson_loglikelihood=0.0;
	int ind2;

	for (unsigned ind =0; ind<contact_value.size();ind++){
		SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);		
		poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distances[ind2]);
		poisson_loglikelihood+=poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter);
	}
	//mexPrintf("contact_sum %f\n",contact_sum);
	poisson_loglikelihood/=contact_sum;
	}
}
 void SIMBA3D::simba3d::compute_poisson_loglikelihood_gradient(){
    if (penalty_weights[0]>0){
	int ind2;
	double scale;
	poisson_loglikelihood=0.0;
	for (unsigned ind =0; ind<contact_value.size();ind++){
		SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);
		poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distances[ind2]);
		poisson_loglikelihood+=poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter);
		scale=-penalty_weights[0]*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter));
		scale/=(distances[ind2]*distances[ind2]);
		scale/=contact_sum;
		x_gradient[contact_row_index[ind]]+=scale*x_differences[ind2];
		y_gradient[contact_row_index[ind]]+=scale*y_differences[ind2];
		z_gradient[contact_row_index[ind]]+=scale*z_differences[ind2];
		x_gradient[contact_column_index[ind]]-=scale*x_differences[ind2];
		y_gradient[contact_column_index[ind]]-=scale*y_differences[ind2];
		z_gradient[contact_column_index[ind]]-=scale*z_differences[ind2];
	}
	poisson_loglikelihood/=contact_sum;
	}
 }

 double SIMBA3D::simba3d::compute_poisson_loglikelihood_gradient_0(double *grad_x,double *grad_y,double *grad_z){
 	poisson_loglikelihood=0.0;
	 if (penalty_weights[0]>0){
		int ind2;

		double scale;

		for (unsigned ind =0; ind<contact_value.size();ind++){
			SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);
			poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distances[ind2]);
			poisson_loglikelihood+=poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter);
			scale=-penalty_weights[0]*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter));
			scale/=(distances[ind2]*distances[ind2]);
			scale/=contact_sum;
			grad_x[contact_row_index[ind]]+=scale*x_differences[ind2];
			grad_y[contact_row_index[ind]]+=scale*y_differences[ind2];
			grad_z[contact_row_index[ind]]+=scale*z_differences[ind2];
			grad_x[contact_column_index[ind]]-=scale*x_differences[ind2];
			grad_y[contact_column_index[ind]]-=scale*y_differences[ind2];
			grad_z[contact_column_index[ind]]-=scale*z_differences[ind2];
		}
	poisson_loglikelihood/=contact_sum;
	}
	return poisson_loglikelihood*penalty_weights[0];
 }
  double SIMBA3D::simba3d::compute_poisson_loglikelihood_gradient_1(double *grad_x,double *grad_y,double *grad_z){
 	poisson_loglikelihood=0.0;
	 if (penalty_weights[0]>0){
		int ind2;

		double scale;
		double x_diff, y_diff, z_diff,distancesqr,distance;
		for (unsigned ind =0; ind<contact_value.size();ind++){
			SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);

			x_diff=x[contact_row_index[ind]]-x[contact_column_index[ind]];
			y_diff=y[contact_row_index[ind]]-y[contact_column_index[ind]];
			z_diff=z[contact_row_index[ind]]-z[contact_column_index[ind]];
			distancesqr=x_diff*x_diff;
			distancesqr+=y_diff*y_diff;
			distancesqr+=z_diff*z_diff;
			distance=std::sqrt(distancesqr);

			poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distance);
			poisson_loglikelihood+=poisson_b_parameter*std::pow(distance,poisson_a_parameter);
			scale=-penalty_weights[0]*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*std::pow(distance,poisson_a_parameter));
			scale/=(distancesqr);
			scale/=contact_sum;
			grad_x[contact_row_index[ind]]+=scale*x_diff;
			grad_y[contact_row_index[ind]]+=scale*y_diff;
			grad_z[contact_row_index[ind]]+=scale*z_diff;
			grad_x[contact_column_index[ind]]-=scale*x_diff;
			grad_y[contact_column_index[ind]]-=scale*y_diff;
			grad_z[contact_column_index[ind]]-=scale*z_diff;
		}
	poisson_loglikelihood/=contact_sum;
	}
	return poisson_loglikelihood*penalty_weights[0];
 }
 void SIMBA3D::simba3d::create_index_mutex(){
	 index_locks.clear();
	 index_locks.resize(30);
	 for (unsigned ii =0 ; ii < index_locks.size();ii++){
		index_locks[ii]=PTHREAD_MUTEX_INITIALIZER;
	 }
 }
 void SIMBA3D::simba3d::compute_poisson_loglikelihood_gradient(double *grad_x,double *grad_y,double *grad_z){
	 if (penalty_weights[0]>0){
		int ind2;
		double scale;
		poisson_loglikelihood=0.0;
		for (unsigned ind =0; ind<contact_value.size();ind++){
			SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);
			poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distances[ind2]);
			poisson_loglikelihood+=poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter);
			scale=-penalty_weights[0]*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter));
			scale/=(distances[ind2]*distances[ind2]);
			scale/=contact_sum;
			//int mutex_id;
			
			//mutex_id=index_locks.size()*(((float) contact_row_index[ind])/ (((float) length)- 1.0));
			//pthread_mutex_lock(&index_locks[mutex_id]);
			grad_x[contact_row_index[ind]]+=scale*x_differences[ind2];
			grad_y[contact_row_index[ind]]+=scale*y_differences[ind2];
			grad_z[contact_row_index[ind]]+=scale*z_differences[ind2];
			//pthread_mutex_unlock(&index_locks[mutex_id]);

			//mutex_id=(((float) contact_column_index[ind])/ ((float) length-1.0));
			//pthread_mutex_lock(&index_locks[mutex_id]);
			grad_x[contact_column_index[ind]]-=scale*x_differences[ind2];
			grad_y[contact_column_index[ind]]-=scale*y_differences[ind2];
			grad_z[contact_column_index[ind]]-=scale*z_differences[ind2];
			//pthread_mutex_unlock(&index_locks[mutex_id]);
		}
		poisson_loglikelihood/=contact_sum;
	 }
 }
 double SIMBA3D::simba3d::compute_poisson_loglikelihood_gradient(int  number_of_threads,double *grad_x,double *grad_y,double *grad_z){
	 if (penalty_weights[0]>0){

		  set_number_of_subtasks(number_of_threads);
		  // initialize threads

		  std::vector<pthread_t> pThreads;
		  pThreads.clear();
		  pThreads.resize(number_of_threads);
		  std::vector<thread_struct> tInfo;
		  tInfo.clear();
		  tInfo.resize(number_of_threads);
		  poisson_loglikelihood=0.0;	
		  for (int ii=0;ii<number_of_threads;ii++){
			tInfo[ii].instance=this;
  			tInfo[ii].subtask_id=ii;
			tInfo[ii].grad_x= grad_x;
			tInfo[ii].grad_y= grad_y;
			tInfo[ii].grad_z= grad_z;
			pthread_create(&pThreads[ii],NULL,split_compute_poisson_loglikelihood_gradient_helper,(void *) &tInfo[ii]);
		  }

		  // wait for all the threads to finish
		  for (int ii=0;ii<number_of_threads;ii++){
			pthread_join(pThreads[ii],NULL);
		  }
		 poisson_loglikelihood/=contact_sum;
		 return poisson_loglikelihood*penalty_weights[3];
	 }
 }
 void SIMBA3D::simba3d::split_compute_poisson_loglikelihood_gradient(int  subtask_id,double *grad_x,double *grad_y,double *grad_z){
	 if (penalty_weights[0]>0){
		int ind2;
		double scale;
			
		for (unsigned ind =subtask_id; ind<contact_value.size();ind+=number_of_subtasks){
			SIMBA3D::triu_ij_to_ind(contact_row_index[ind],contact_column_index[ind],ind2);
			poisson_loglikelihood-=poisson_a_parameter*contact_value[ind]*log(distances[ind2]);
			poisson_loglikelihood+=poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter);
			scale=-penalty_weights[0]*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*std::pow(distances[ind2],poisson_a_parameter));
			scale/=(distances[ind2]*distances[ind2]);
			scale/=contact_sum;
			int mutex_id;
			
			mutex_id=index_locks.size()*(((float) contact_row_index[ind])/ (((float) length)- 1.0));
			pthread_mutex_lock(&index_locks[mutex_id]);
			grad_x[contact_row_index[ind]]+=scale*x_differences[ind2];
			grad_y[contact_row_index[ind]]+=scale*y_differences[ind2];
			grad_z[contact_row_index[ind]]+=scale*z_differences[ind2];
			pthread_mutex_unlock(&index_locks[mutex_id]);

			mutex_id=(((float) contact_column_index[ind])/ ((float) length-1.0));
			pthread_mutex_lock(&index_locks[mutex_id]);
			grad_x[contact_column_index[ind]]-=scale*x_differences[ind2];
			grad_y[contact_column_index[ind]]-=scale*y_differences[ind2];
			grad_z[contact_column_index[ind]]-=scale*z_differences[ind2];
			pthread_mutex_unlock(&index_locks[mutex_id]);
		}
		
	 }
 }
//.................................................................................
void SIMBA3D::simba3d::compute_h1_gradient()
{	
	if (penalty_weights[1]>0){
	double L_squared=sum_adjacent_distance*sum_adjacent_distance;
	double L_cubed=sum_adjacent_distance*sum_adjacent_distance*sum_adjacent_distance;
	double S=sum_of_squared_adjacent_distance;
	double term0=2.0/L_squared;
	double coef=(length-1)*(term0-2.0*S/(L_cubed*distances[0]));
	x_gradient[0]+=penalty_weights[1]*coef*x_differences[0];
	y_gradient[0]+=penalty_weights[1]*coef*y_differences[0];
	z_gradient[0]+=penalty_weights[1]*coef*z_differences[0];
	int adj_ind,adj_ind0,ii;
	double coef0;
	for (ii=1;ii<length-1;ii++){
		adj_ind=(ii+2)*(ii+1)/2 -1;// get the adjacent index
		adj_ind0=(ii+1)*(ii)/2 -1;// get the previous adjacent index
		coef0=coef;
		//coef=(length-1)*(term0+2.0*S/(L_cubed*distances[adj_ind0]));
		x_gradient[ii]-=penalty_weights[1]*x_differences[adj_ind0]*coef0;
		y_gradient[ii]-=penalty_weights[1]*y_differences[adj_ind0]*coef0;
		z_gradient[ii]-=penalty_weights[1]*z_differences[adj_ind0]*coef0;
		coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
		x_gradient[ii]+=penalty_weights[1]*x_differences[adj_ind]*coef;
		y_gradient[ii]+=penalty_weights[1]*y_differences[adj_ind]*coef;
		z_gradient[ii]+=penalty_weights[1]*z_differences[adj_ind]*coef;
	}
	ii=length-2;
	adj_ind=(ii+2)*(ii+1)/2 -1;// get the adjacent index
	coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
	x_gradient[length-1]-=penalty_weights[1]*x_differences[adj_ind]*coef;
	y_gradient[length-1]-=penalty_weights[1]*y_differences[adj_ind]*coef;
	z_gradient[length-1]-=penalty_weights[1]*z_differences[adj_ind]*coef;
	}
}
double SIMBA3D::simba3d::compute_h1_gradient(double *grad_x,double *grad_y,double *grad_z)
{	if (penalty_weights[1]>0){
	double L_squared=sum_adjacent_distance*sum_adjacent_distance;
	double L_cubed=sum_adjacent_distance*sum_adjacent_distance*sum_adjacent_distance;
	double S=sum_of_squared_adjacent_distance;
	double term0=2.0/L_squared;
	double coef=(length-1)*(term0-2.0*S/(L_cubed*distances[0]));
	grad_x[0]+=penalty_weights[1]*coef*x_differences[0];
	grad_y[0]+=penalty_weights[1]*coef*y_differences[0];
	grad_z[0]+=penalty_weights[1]*coef*z_differences[0];
	int adj_ind,adj_ind0,ii;
	double coef0;
	for (ii=1;ii<length-1;ii++){
		adj_ind=(ii+2)*(ii+1)/2 -1;// get the adjacent index
		adj_ind0=(ii+1)*(ii)/2 -1;// get the previous adjacent index
		coef0=coef;
		//coef=(length-1)*(term0+2.0*S/(L_cubed*distances[adj_ind0]));
		grad_x[ii]-=penalty_weights[1]*x_differences[adj_ind0]*coef0;
		grad_y[ii]-=penalty_weights[1]*y_differences[adj_ind0]*coef0;
		grad_z[ii]-=penalty_weights[1]*z_differences[adj_ind0]*coef0;
		coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
		grad_x[ii]+=penalty_weights[1]*x_differences[adj_ind]*coef;
		grad_y[ii]+=penalty_weights[1]*y_differences[adj_ind]*coef;
		grad_z[ii]+=penalty_weights[1]*z_differences[adj_ind]*coef;
	}
	ii=length-2;
	adj_ind=(ii+2)*(ii+1)/2 -1;// get the adjacent index
	coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
	grad_x[length-1]-=penalty_weights[1]*x_differences[adj_ind]*coef;
	grad_y[length-1]-=penalty_weights[1]*y_differences[adj_ind]*coef;
	grad_z[length-1]-=penalty_weights[1]*z_differences[adj_ind]*coef;
}
	return penalty_weights[1]*h1_uniform_spacing_value;
}
double SIMBA3D::simba3d::compute_h2a_gradient(double *grad_x,double *grad_y,double *grad_z)
{
	
	if (penalty_weights[4]>0){
		h2a_value=0.0;
		int adj_ind0,adj_ind1,adj_ind2,adj_ind3,adj_indn1,adj_indn2,adj_indn3,adj_indn4;
		double dot,coef,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9;
		///
		int ii;
		for (ii=1; ii<length-1;ii++){
			adj_ind0=floor(((ii+1)*(ii))/2.0 -1);// get the adjacent indexes
			adj_ind1=floor(((ii+2)*(ii+1))/2.0 -1);// get the adjacent indexes
			dot =(x_differences[adj_ind0]*x_differences[adj_ind1]);
			dot +=(y_differences[adj_ind0]*y_differences[adj_ind1]);
			dot +=(z_differences[adj_ind0]*z_differences[adj_ind1]);
			h2a_value-=dot/(distances[adj_ind0]*distances[adj_ind1]);	
			
			if (ii==1){
				// compute the first gradient term
				coef0 =-x_differences[adj_ind0]*x_differences[adj_ind1];
				coef0+=-y_differences[adj_ind0]*y_differences[adj_ind1];
				coef0+=-z_differences[adj_ind0]*z_differences[adj_ind1];
				coef0/=pow(distances[adj_ind0],3.0)*distances[adj_ind1]	;	

				grad_x[0]-=(penalty_weights[4]*x_differences[adj_ind0]*coef0)/(length-2);
				grad_y[0]-=(penalty_weights[4]*y_differences[adj_ind0]*coef0)/(length-2);
				grad_z[0]-=(penalty_weights[4]*z_differences[adj_ind0]*coef0)/(length-2);
				/*
				std::cout <<coef0<<"cpp\n";
				std::cout <<x_differences[adj_ind0]*coef0<<"cpp\n";
				std::cout <<y_differences[adj_ind0]*coef0<<"cpp\n";
				std::cout <<z_differences[adj_ind0]*coef0<<"cpp\n";
				*/
				coef=(distances[adj_ind0]*distances[adj_ind1]);

				grad_x[0]-=(penalty_weights[4]*x_differences[adj_ind1]/coef)/(length-2);
				grad_y[0]-=(penalty_weights[4]*y_differences[adj_ind1]/coef)/(length-2);
				grad_z[0]-=(penalty_weights[4]*z_differences[adj_ind1]/coef)/(length-2);		
				//std::cout << grad_x[0] <<" " << grad_y[0] <<" "<< grad_z[0] <<"cpp\n ";

			} else if (ii==2) {
				// compute the second gradient term
				adj_ind0=floor(((2)*(1))/2.0 -1);// get the adjacent indexes;
				adj_ind1=floor(((3)*(2))/2.0 -1);// get the adjacent indexes;
				adj_ind2=floor(((4)*(3))/2.0 -1);// get the adjacent indexes;
				coef0 =x_differences[adj_ind0]*x_differences[adj_ind1];
				coef0+=y_differences[adj_ind0]*y_differences[adj_ind1];
				coef0+=z_differences[adj_ind0]*z_differences[adj_ind1];
				coef1 =(pow(distances[adj_ind0],3.0)*distances[adj_ind1]);	
				coef2 =1/(distances[adj_ind0]*distances[adj_ind1]);
				coef  =(coef0/coef1+coef2);
				grad_x[1]-=penalty_weights[4]*x_differences[adj_ind0]*coef/(length-2);
				grad_y[1]-=penalty_weights[4]*y_differences[adj_ind0]*coef/(length-2);
				grad_z[1]-=penalty_weights[4]*z_differences[adj_ind0]*coef/(length-2);

				coef3 =(distances[adj_ind0]*pow(distances[adj_ind1],3))	;
				coef4 =x_differences[adj_ind1]*x_differences[adj_ind2];
				coef4+=y_differences[adj_ind1]*y_differences[adj_ind2];
				coef4+=z_differences[adj_ind1]*z_differences[adj_ind2];
				coef5 =(pow(distances[adj_ind1],3.0)*distances[adj_ind2]);
				coef  =(-coef2-coef0/coef3- coef4/coef5);
				grad_x[1]+=-penalty_weights[4]*x_differences[adj_ind1]*coef/(length-2);
				grad_y[1]+=-penalty_weights[4]*y_differences[adj_ind1]*coef/(length-2);
				grad_z[1]+=-penalty_weights[4]*z_differences[adj_ind1]*coef/(length-2);

				coef6  =1/(distances[adj_ind1]*distances[adj_ind2]);
				grad_x[1]+=-penalty_weights[4]*x_differences[adj_ind2]*coef6/(length-2);
				grad_y[1]+=-penalty_weights[4]*y_differences[adj_ind2]*coef6/(length-2);
				grad_z[1]+=-penalty_weights[4]*z_differences[adj_ind2]*coef6/(length-2);

			} else{
				adj_indn3=floor(((ii+-1)*(ii-2))/2.0 -1);// get the adjacent indexes
				adj_indn2=floor(((ii+0)*(ii-1))/2.0 -1);// get the adjacent indexes
				adj_indn1=floor(((ii+1)*(ii+0))/2.0 -1);// get the adjacent indexes
				adj_ind0=floor(((ii+2)*(ii+1))/2.0 -1);// get the adjacent indexes
				coef0  =1/(distances[adj_indn3]*distances[adj_indn2]);
				grad_x[ii-1]+=penalty_weights[4]*x_differences[adj_indn3]*coef0/(length-2);
				grad_y[ii-1]+=penalty_weights[4]*y_differences[adj_indn3]*coef0/(length-2);
				grad_z[ii-1]+=penalty_weights[4]*z_differences[adj_indn3]*coef0/(length-2);
				
				coef1 =x_differences[adj_indn3]*x_differences[adj_indn2];
				coef1+=y_differences[adj_indn3]*y_differences[adj_indn2];
				coef1+=z_differences[adj_indn3]*z_differences[adj_indn2];
				coef2 =(distances[adj_indn3]*pow(distances[adj_indn2],3.0));
				coef3 =1/(distances[adj_indn2]*distances[adj_indn1]);
				coef4 =x_differences[adj_indn2]*x_differences[adj_indn1];
				coef4+=y_differences[adj_indn2]*y_differences[adj_indn1];
				coef4+=z_differences[adj_indn2]*z_differences[adj_indn1];
				coef5=pow(distances[adj_indn2],3.0)*distances[adj_indn1];
				coef=-(coef1/coef2+coef3+coef4/coef5);
				grad_x[ii-1]+=penalty_weights[4]*x_differences[adj_indn2]*coef/(length-2);
				grad_y[ii-1]+=penalty_weights[4]*y_differences[adj_indn2]*coef/(length-2);
				grad_z[ii-1]+=penalty_weights[4]*z_differences[adj_indn2]*coef/(length-2);	

					
				coef6=(distances[adj_indn2]*pow(distances[adj_indn1],3.0))	;	
				coef7 =x_differences[adj_indn1]*x_differences[adj_ind0];
				coef7+=y_differences[adj_indn1]*y_differences[adj_ind0];
				coef7+=z_differences[adj_indn1]*z_differences[adj_ind0]	;	
				coef8=pow(distances[adj_indn1],3.0)*distances[adj_ind0];
				coef9=pow(distances[adj_indn1],3.0)*distances[adj_indn2];			
				coef=coef3 +coef4/coef9+coef7/coef8;
				grad_x[ii-1]+=penalty_weights[4]*x_differences[adj_indn1]*coef/(length-2);
				grad_y[ii-1]+=penalty_weights[4]*y_differences[adj_indn1]*coef/(length-2);
				grad_z[ii-1]+=penalty_weights[4]*z_differences[adj_indn1]*coef/(length-2);

				coef=-1/(distances[adj_indn1]*distances[adj_ind0]);
				grad_x[ii-1]+=penalty_weights[4]*x_differences[adj_ind0]*coef/(length-2);
				grad_y[ii-1]+=penalty_weights[4]*y_differences[adj_ind0]*coef/(length-2);
				grad_z[ii-1]+=penalty_weights[4]*z_differences[adj_ind0]*coef/(length-2);	
							
			}
		}	
		/*
		ii=3;	

		adj_indn3=floor(((ii+-1)*(ii-2))/2.0 -1);// get the adjacent indexes
		adj_indn2=floor(((ii+0)*(ii-1))/2.0 -1);// get the adjacent indexes
		adj_indn1=floor(((ii+1)*(ii+0))/2.0 -1);// get the adjacent indexes
		adj_ind0=floor(((ii+2)*(ii+1))/2.0 -1);// get the adjacent indexes
		coef0  =1/(distances[adj_indn3]*distances[adj_indn2]);		
		coef1 =x_differences[adj_indn3]*x_differences[adj_indn2];
		coef1+=y_differences[adj_indn3]*y_differences[adj_indn2];
		coef1+=z_differences[adj_indn3]*z_differences[adj_indn2];
		coef2 =(distances[adj_indn3]*pow(distances[adj_indn2],3.0));
		coef3 =1/(distances[adj_indn2]*distances[adj_indn1]);
		coef4 =x_differences[adj_indn2]*x_differences[adj_indn1];
		coef4+=y_differences[adj_indn2]*y_differences[adj_indn1];
		coef4+=z_differences[adj_indn2]*z_differences[adj_indn1];
		coef5=pow(distances[adj_indn2],3.0)*distances[adj_indn1];
		coef=-(coef1/coef2+coef3+coef4/coef5);
		coef6=(distances[adj_indn2]*pow(distances[adj_indn1],3.0))	;	
		coef7 =x_differences[adj_indn1]*x_differences[adj_ind0];
		coef7+=y_differences[adj_indn1]*y_differences[adj_ind0];
		coef7+=z_differences[adj_indn1]*z_differences[adj_ind0]	;	
		coef8=pow(distances[adj_indn1],3.0)*distances[adj_ind0];
		coef9=pow(distances[adj_indn1],3.0)*distances[adj_indn2];
		coef=coef3 +coef4/coef9+coef7/coef8;
		std::cout << x_differences[adj_indn1]*coef<<" "<<y_differences[adj_indn1]*coef<<" " << z_differences[adj_indn1]*coef <<"cpp\n ";
		*/
		

		//
		// compute the next to last gradient 	term
		adj_indn4=floor(((length+-2)*(length-3))/2.0 -1);// get the adjacent indexes
		adj_indn3=floor(((length+-1)*(length-2))/2.0 -1);// get the adjacent indexes
		adj_indn2=floor(((length+-0)*(length-1))/2.0 -1);// get the adjacent indexes	
		coef=1/(distances[adj_indn4]*distances[adj_indn3]);

		grad_x[length-2]+=penalty_weights[4]*x_differences[adj_indn4]*coef/(length-2);
		grad_y[length-2]+=penalty_weights[4]*y_differences[adj_indn4]*coef/(length-2);
		grad_z[length-2]+=penalty_weights[4]*z_differences[adj_indn4]*coef/(length-2);
		//
		coef1 =x_differences[adj_indn4]*x_differences[adj_indn3];
		coef1+=y_differences[adj_indn4]*y_differences[adj_indn3];
		coef1+=z_differences[adj_indn4]*z_differences[adj_indn3];
		coef2=pow(distances[adj_indn3],3.0)*distances[adj_indn4];
		coef3=1/(distances[adj_indn3]*distances[adj_indn2]);
		coef4 =x_differences[adj_indn3]*x_differences[adj_indn2];
		coef4+=y_differences[adj_indn3]*y_differences[adj_indn2];
		coef4+=z_differences[adj_indn3]*z_differences[adj_indn2];	
		coef5=(pow(distances[adj_indn3],3.0)*distances[adj_indn2]);
		coef=-(coef1/coef2+coef3+coef4/coef5);

		grad_x[length-2]+=penalty_weights[4]*x_differences[adj_indn3]*coef/(length-2);
		grad_y[length-2]+=penalty_weights[4]*y_differences[adj_indn3]*coef/(length-2);
		grad_z[length-2]+=penalty_weights[4]*z_differences[adj_indn3]*coef/(length-2);
		//
		coef6=(distances[adj_indn3]*pow(distances[adj_indn2],3.0));
		coef=(-coef3- coef4/coef6);
		//coef=(-coef3- coef4/(w[n-3]*pow(w[n-2],3.0)))
		grad_x[length-2]-=penalty_weights[4]*x_differences[adj_indn2]*coef/(length-2);
		grad_y[length-2]-=penalty_weights[4]*y_differences[adj_indn2]*coef/(length-2);
		grad_z[length-2]-=penalty_weights[4]*z_differences[adj_indn2]*coef/(length-2);	

		grad_x[length-1]+=penalty_weights[4]*x_differences[adj_indn3]*coef3/(length-2);
		grad_y[length-1]+=penalty_weights[4]*y_differences[adj_indn3]*coef3/(length-2);
		grad_z[length-1]+=penalty_weights[4]*z_differences[adj_indn3]*coef3/(length-2);
		/*
		std::cout << x_differences[adj_indn3]*coef3 <<" "<<"cpp\n ";
		std::cout << y_differences[adj_indn3]*coef3 <<" "<<"cpp\n ";
		std::cout << z_differences[adj_indn3]*coef3 <<" "<<"cpp\n ";
		*/
		coef=(coef4/coef6);
		grad_x[length-1]-=penalty_weights[4]*x_differences[adj_indn2]*coef/(length-2);
		grad_y[length-1]-=penalty_weights[4]*y_differences[adj_indn2]*coef/(length-2);
		grad_z[length-1]-=penalty_weights[4]*z_differences[adj_indn2]*coef/(length-2);			
		//////////////////
		/*	
		std::cout << -x_differences[adj_indn2]*coef <<" "<<"cpp\n ";
		std::cout << -y_differences[adj_indn2]*coef <<" "<<"cpp\n ";
		std::cout << -z_differences[adj_indn2]*coef <<" "<<"cpp\n ";
		*/

		h2a_value/=length-2	;
	/*


	*/
	} 

	return h2a_value*penalty_weights[4];
}
double SIMBA3D::simba3d::compute_h2a_gradient()
{
	return compute_h2a_gradient(&x_gradient[0],&y_gradient[0],&z_gradient[0]);
}
//.................................................................................
void SIMBA3D::simba3d::compute_h2b_computations()
{
	if ((penalty_weights[2]>0)|(penalty_weights[3]>0)){
	sum_pairwise_distance=0.0;
	if (F_h2b.size()!= (unsigned)number_of_pairs){
		F_h2b.clear();
		F_h2b.reserve(number_of_pairs);
		F_h2b.resize(number_of_pairs);
	}
	if (H_h2b.size()!=(unsigned) number_of_pairs){
		H_h2b.clear();
		H_h2b.reserve(number_of_pairs);
		H_h2b.resize(number_of_pairs);
	}
	if (G_h2b.size()!= (unsigned) number_of_pairs){
		G_h2b.clear();
		G_h2b.reserve(number_of_pairs);
		G_h2b.resize(number_of_pairs);
	}	
	h2b_pairwise_spacing_value=0.0;
	int ind;
	for (ind=0;ind<number_of_pairs;ind++){
		//SIMBA3D::triu_ind_to_ij(ii ,jj,ind);
		double u_jk= distances[ind]/average_adjacent_distance;
		F_h2b[ind]=1.0-u_jk;

		G_h2b[ind]=exp(alpha_h2b*(u_jk-1.0));
		H_h2b[ind]=1.0/(1.0+G_h2b[ind]);
		h2b_pairwise_spacing_value+=F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind];
	}
	h2b_pairwise_spacing_value*=2.0/((double) (length*(length-1.0)));
	}
}

//.................................................................................
void SIMBA3D::simba3d::compute_h2b_gradient()
{
	  if (penalty_weights[2]>0){
	int ind,ind2,jj,kk;
	for (int ii=0; ii <length; ii++){
		// partial derivative of mu
		double partial_mu_xi_x;
		double partial_mu_xi_y;
		double partial_mu_xi_z;
		if (ii==0){
			SIMBA3D::triu_ij_to_ind(0,1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
		} else if (ii==length-1){
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x=-x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=-y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=-z_differences[ind2]/((length-1)*distances[ind2]);
		} else{
			SIMBA3D::triu_ij_to_ind(ii,ii+1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x-=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y-=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z-=z_differences[ind2]/((length-1)*distances[ind2]);
		}		
		for (ind=0;ind<number_of_pairs;ind++){
			SIMBA3D::triu_ind_to_ij(jj,kk,ind);
			// partial derivative of distances
			double partial_djk_xi_x=0.0;
			double partial_djk_xi_y=0.0;
			double partial_djk_xi_z=0.0;
			if (jj==ii){
				partial_djk_xi_x=x_differences[ind]/distances[ind];
				partial_djk_xi_y=y_differences[ind]/distances[ind];
				partial_djk_xi_z=z_differences[ind]/distances[ind];
			} else if (kk==ii){
				partial_djk_xi_x=-x_differences[ind]/distances[ind];
				partial_djk_xi_y=-y_differences[ind]/distances[ind];
				partial_djk_xi_z=-z_differences[ind]/distances[ind];
			}
			// partial derivative of ujk
			double partial_ujk_xi_x=(partial_djk_xi_x+(F_h2b[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
			double partial_ujk_xi_y=(partial_djk_xi_y+(F_h2b[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
			double partial_ujk_xi_z=(partial_djk_xi_z+(F_h2b[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
			// compute gradient
			double scale=-penalty_weights[2]*2/(length*(length-1));
			double term=scale*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind]*(4.0+alpha_h2b*F_h2b[ind]*G_h2b[ind]*H_h2b[ind]);
			x_gradient[ii]+=term*partial_ujk_xi_x;
			y_gradient[ii]+=term*partial_ujk_xi_y;
			z_gradient[ii]+=term*partial_ujk_xi_z;
			/*
			// partial derivative of Fjk
			double partial_Fjk_xi_x=-partial_ujk_xi_x;
			double partial_Fjk_xi_y=-partial_ujk_xi_y;
			double partial_Fjk_xi_z=-partial_ujk_xi_z;
			// partial derivative of Hjk
			double partial_Hjk_xi_x=-alpha_h2b*partial_ujk_xi_x*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_y=-alpha_h2b*partial_ujk_xi_y*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_z=-alpha_h2b*partial_ujk_xi_z*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			// compute gradient
			double term1=scale*4.0*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind];
			double term2=scale*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*F_h2b[ind];
			x_gradient[ii]+=term1*partial_Fjk_xi_x+term2*partial_Hjk_xi_x;
			y_gradient[ii]+=term1*partial_Fjk_xi_y+term2*partial_Hjk_xi_y;
			z_gradient[ii]+=term1*partial_Fjk_xi_z+term2*partial_Hjk_xi_z;
			*/

			
		}
	}
	  }
}
double SIMBA3D::simba3d::compute_h2b_gradient(double *grad_x,double *grad_y,double *grad_z)
{
	if (penalty_weights[2]>0){
		int ind,ind2,jj,kk;
		for (int ii=0; ii <length; ii++){
			// partial derivative of mu
			double partial_mu_xi_x;
			double partial_mu_xi_y;
			double partial_mu_xi_z;
			if (ii==0){
				SIMBA3D::triu_ij_to_ind(0,1,ind2);
				partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
			} else if (ii==length-1){
				SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
				partial_mu_xi_x=-x_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_y=-y_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_z=-z_differences[ind2]/((length-1)*distances[ind2]);
			} else{
				SIMBA3D::triu_ij_to_ind(ii,ii+1,ind2);
				partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
				SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
				partial_mu_xi_x-=x_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_y-=y_differences[ind2]/((length-1)*distances[ind2]);
				partial_mu_xi_z-=z_differences[ind2]/((length-1)*distances[ind2]);
			}		
			for (ind=0;ind<number_of_pairs;ind++){
				SIMBA3D::triu_ind_to_ij(jj,kk,ind);
				// partial derivative of distances
				double partial_djk_xi_x=0.0;
				double partial_djk_xi_y=0.0;
				double partial_djk_xi_z=0.0;
				if (jj==ii){
					partial_djk_xi_x=x_differences[ind]/distances[ind];
					partial_djk_xi_y=y_differences[ind]/distances[ind];
					partial_djk_xi_z=z_differences[ind]/distances[ind];
				} else if (kk==ii){
					partial_djk_xi_x=-x_differences[ind]/distances[ind];
					partial_djk_xi_y=-y_differences[ind]/distances[ind];
					partial_djk_xi_z=-z_differences[ind]/distances[ind];
				}
				// partial derivative of ujk
				double partial_ujk_xi_x=(partial_djk_xi_x+(F_h2b[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
				double partial_ujk_xi_y=(partial_djk_xi_y+(F_h2b[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
				double partial_ujk_xi_z=(partial_djk_xi_z+(F_h2b[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
				// compute gradient
				double scale=-penalty_weights[2]*2/(length*(length-1));
				double term=scale*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind]*(4.0+alpha_h2b*F_h2b[ind]*G_h2b[ind]*H_h2b[ind]);
				grad_x[ii]+=term*partial_ujk_xi_x;
				grad_y[ii]+=term*partial_ujk_xi_y;
				grad_z[ii]+=term*partial_ujk_xi_z;
			}
		}
	}
	return h2b_pairwise_spacing_value*penalty_weights[2];
}
//.................................................................................
void SIMBA3D::simba3d::compute_h2c_computations(){
	if (penalty_weights[3]>0){
	h2c_pairwise_spacing_value=0.0;
	int ind;	
	for (ind=0;ind<number_of_pairs;ind++)
	{
		if (distances[ind]/average_adjacent_distance<h2c_radius){
			h2c_pairwise_spacing_value+=std::pow(h2c_radius-distances[ind]/average_adjacent_distance,4.0);
		}
	}	
	h2c_pairwise_spacing_value*=2.0/((double) (length*(length-1.0)));
	}
}
void SIMBA3D::simba3d::compute_h2c_gradient(){
if (penalty_weights[3]>0){
	int ind,ind2,jj,kk;
	for (int ii=0; ii <length; ii++){
		// partial derivative of mu
		double partial_mu_xi_x;
		double partial_mu_xi_y;
		double partial_mu_xi_z;
		if (ii==0){
			SIMBA3D::triu_ij_to_ind(0,1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
		} else if (ii==length-1){
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x=-x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=-y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=-z_differences[ind2]/((length-1)*distances[ind2]);
		} else{
			SIMBA3D::triu_ij_to_ind(ii,ii+1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x-=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y-=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z-=z_differences[ind2]/((length-1)*distances[ind2]);
		}		
		for (ind=0;ind<number_of_pairs;ind++){

			SIMBA3D::triu_ind_to_ij(jj,kk,ind);
			if ((1.0-F_h2b[ind])<h2c_radius){
				// partial derivative of distances
				double partial_djk_xi_x=0.0;
				double partial_djk_xi_y=0.0;
				double partial_djk_xi_z=0.0;
				if (jj==ii){
					partial_djk_xi_x=x_differences[ind]/distances[ind];
					partial_djk_xi_y=y_differences[ind]/distances[ind];
					partial_djk_xi_z=z_differences[ind]/distances[ind];
				} else if (kk==ii){
					partial_djk_xi_x=-x_differences[ind]/distances[ind];
					partial_djk_xi_y=-y_differences[ind]/distances[ind];
					partial_djk_xi_z=-z_differences[ind]/distances[ind];
				}
				// partial derivative of ujk
				double partial_ujk_xi_x=(partial_djk_xi_x+(F_h2b[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
				double partial_ujk_xi_y=(partial_djk_xi_y+(F_h2b[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
				double partial_ujk_xi_z=(partial_djk_xi_z+(F_h2b[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
				// compute gradient
				double scale=-penalty_weights[3]*2/(length*(length-1));
				double F=(h2c_radius-distances[ind]/average_adjacent_distance);
				double term=scale*4*F*F*F;
				
				x_gradient[ii]+=term*partial_ujk_xi_x;
				y_gradient[ii]+=term*partial_ujk_xi_y;
				z_gradient[ii]+=term*partial_ujk_xi_z;
			}
			/*
			// partial derivative of Fjk
			double partial_Fjk_xi_x=-partial_ujk_xi_x;
			double partial_Fjk_xi_y=-partial_ujk_xi_y;
			double partial_Fjk_xi_z=-partial_ujk_xi_z;
			// partial derivative of Hjk
			double partial_Hjk_xi_x=-alpha_h2b*partial_ujk_xi_x*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_y=-alpha_h2b*partial_ujk_xi_y*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_z=-alpha_h2b*partial_ujk_xi_z*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			// compute gradient
			double term1=scale*4.0*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind];
			double term2=scale*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*F_h2b[ind];
			x_gradient[ii]+=term1*partial_Fjk_xi_x+term2*partial_Hjk_xi_x;
			y_gradient[ii]+=term1*partial_Fjk_xi_y+term2*partial_Hjk_xi_y;
			z_gradient[ii]+=term1*partial_Fjk_xi_z+term2*partial_Hjk_xi_z;
			*/

			
		}
	}
}
}
void SIMBA3D::simba3d::compute_h2c_gradient(double *grad_x,double *grad_y,double *grad_z){
	if (penalty_weights[3]>0){
	int ind,ind2,jj,kk;
	double scale=-penalty_weights[3]*4.0*(2.0)/((double) length*(length-1));
	for (int ii=0; ii <length; ii++){
		// partial derivative of mu
		double partial_mu_xi_x;
		double partial_mu_xi_y;
		double partial_mu_xi_z;
		if (ii==0){
			SIMBA3D::triu_ij_to_ind(0,1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
		} else if (ii==length-1){
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x=-x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=-y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=-z_differences[ind2]/((length-1)*distances[ind2]);
		} else{
			SIMBA3D::triu_ij_to_ind(ii,ii+1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x-=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y-=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z-=z_differences[ind2]/((length-1)*distances[ind2]);
		}	

		
		for (ind=0;ind<number_of_pairs;ind++){

			SIMBA3D::triu_ind_to_ij(jj,kk,ind);
			if ((1.0-F_h2b[ind])<h2c_radius){
				// partial derivative of distances
				double partial_djk_xi_x=0.0;
				double partial_djk_xi_y=0.0;
				double partial_djk_xi_z=0.0;
				if (jj==ii){
					partial_djk_xi_x=x_differences[ind]/distances[ind];
					partial_djk_xi_y=y_differences[ind]/distances[ind];
					partial_djk_xi_z=z_differences[ind]/distances[ind];
				} else if (kk==ii){
					partial_djk_xi_x=-x_differences[ind]/distances[ind];
					partial_djk_xi_y=-y_differences[ind]/distances[ind];
					partial_djk_xi_z=-z_differences[ind]/distances[ind];
				}
				// partial derivative of ujk
				double partial_ujk_xi_x=(partial_djk_xi_x+(F_h2b[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
				double partial_ujk_xi_y=(partial_djk_xi_y+(F_h2b[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
				double partial_ujk_xi_z=(partial_djk_xi_z+(F_h2b[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
				//mexPrintf("%.16f\t",partial_ujk_xi_x);
				//mexPrintf("%.16f\t",partial_ujk_xi_y);
				// compute gradient
				
				double F=(h2c_radius-distances[ind]/average_adjacent_distance);
				double term=scale*exp(3.0*log(F));
				
				grad_x[ii]+=term*partial_ujk_xi_x;
				grad_y[ii]+=term*partial_ujk_xi_y;
				grad_z[ii]+=term*partial_ujk_xi_z;
			}
			/*
			// partial derivative of Fjk
			double partial_Fjk_xi_x=-partial_ujk_xi_x;
			double partial_Fjk_xi_y=-partial_ujk_xi_y;
			double partial_Fjk_xi_z=-partial_ujk_xi_z;
			// partial derivative of Hjk
			double partial_Hjk_xi_x=-alpha_h2b*partial_ujk_xi_x*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_y=-alpha_h2b*partial_ujk_xi_y*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			double partial_Hjk_xi_z=-alpha_h2b*partial_ujk_xi_z*G_h2b[ind]*H_h2b[ind]*H_h2b[ind];
			// compute gradient
			double term1=scale*4.0*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*H_h2b[ind];
			double term2=scale*F_h2b[ind]*F_h2b[ind]*F_h2b[ind]*F_h2b[ind];
			x_gradient[ii]+=term1*partial_Fjk_xi_x+term2*partial_Hjk_xi_x;
			y_gradient[ii]+=term1*partial_Fjk_xi_y+term2*partial_Hjk_xi_y;
			z_gradient[ii]+=term1*partial_Fjk_xi_z+term2*partial_Hjk_xi_z;
			*/

			
		}
	}
	}
}
double SIMBA3D::simba3d::compute_h2c_gradient(int number_of_threads,double *grad_x,double *grad_y,double *grad_z)
{
  set_number_of_subtasks(number_of_threads);
  // initialize threads

  std::vector<pthread_t> pThreads;
  pThreads.clear();
  pThreads.resize(number_of_threads);
  std::vector<thread_struct> tInfo;
  tInfo.clear();
  tInfo.resize(number_of_threads);
  for (int ii=0;ii<number_of_threads;ii++){

	tInfo[ii].instance=this;
	tInfo[ii].subtask_id=ii;
	tInfo[ii].grad_x= grad_x;
	tInfo[ii].grad_y= grad_y;
	tInfo[ii].grad_z= grad_z;

	pthread_create(&pThreads[ii],NULL,split_compute_h2c_gradient_helper,(void *) &tInfo[ii]);
  }
  // wait for all the threads to finish
  for (int ii=0;ii<number_of_threads;ii++){
	pthread_join(pThreads[ii],NULL);
  }
  return h2c_pairwise_spacing_value*penalty_weights[3];
}
void SIMBA3D::simba3d::split_compute_h2c_gradient(int subtask_id,double *grad_x,double *grad_y,double *grad_z)
{
	int ind,ind2,jj,kk;
	//std::cout <<length<<"started\n";
	double scale=-penalty_weights[3]*4.0*(2.0)/((double) length*(length-1));
	for (int ii=subtask_id; ii <length; ii+=number_of_subtasks){
		// partial derivative of mu
		  	  	

		double partial_mu_xi_x;
		double partial_mu_xi_y;
		double partial_mu_xi_z;
		if (ii==0){
			SIMBA3D::triu_ij_to_ind(0,1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
		} else if (ii==length-1){
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x=-x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=-y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=-z_differences[ind2]/((length-1)*distances[ind2]);
		} else{
			SIMBA3D::triu_ij_to_ind(ii,ii+1,ind2);
			partial_mu_xi_x=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z=z_differences[ind2]/((length-1)*distances[ind2]);
			SIMBA3D::triu_ij_to_ind(ii-1,ii,ind2);
			partial_mu_xi_x-=x_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_y-=y_differences[ind2]/((length-1)*distances[ind2]);
			partial_mu_xi_z-=z_differences[ind2]/((length-1)*distances[ind2]);
		}			
		for (ind=0;ind<number_of_pairs;ind++){

			SIMBA3D::triu_ind_to_ij(jj,kk,ind);
			if ((1.0-F_h2b[ind])<h2c_radius){
				// partial derivative of distances
				double partial_djk_xi_x=0.0;
				double partial_djk_xi_y=0.0;
				double partial_djk_xi_z=0.0;
				if (jj==ii){
					partial_djk_xi_x=x_differences[ind]/distances[ind];
					partial_djk_xi_y=y_differences[ind]/distances[ind];
					partial_djk_xi_z=z_differences[ind]/distances[ind];
				} else if (kk==ii){
					partial_djk_xi_x=-x_differences[ind]/distances[ind];
					partial_djk_xi_y=-y_differences[ind]/distances[ind];
					partial_djk_xi_z=-z_differences[ind]/distances[ind];
				}
				// partial derivative of ujk
				double partial_ujk_xi_x=(partial_djk_xi_x+(F_h2b[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
				double partial_ujk_xi_y=(partial_djk_xi_y+(F_h2b[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
				double partial_ujk_xi_z=(partial_djk_xi_z+(F_h2b[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
				//mexPrintf("%.16f\t",partial_ujk_xi_x);
				//mexPrintf("%.16f\t",partial_ujk_xi_y);
				// compute gradient
				
				double F=(h2c_radius-distances[ind]/average_adjacent_distance);
				double term=scale*exp(3.0*log(F));
				
				grad_x[ii]+=term*partial_ujk_xi_x;
				grad_y[ii]+=term*partial_ujk_xi_y;
				grad_z[ii]+=term*partial_ujk_xi_z;
			}

			
		}
	}

}
double SIMBA3D::simba3d::compute_energy(){
  //compute_pairwise_computations();
  //compute_adjacent_computations();
  
  //compute_poisson_loglikelihood();
  /*
  if (penalty_weights[0]>0){
	  start = std::clock();
      compute_poisson_loglikelihood();
	  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	  mexPrintf("\tcompute_poisson_loglikelihood %f\n",duration);
  }
  if ((penalty_weights[2]>0)|(penalty_weights[3]>0)){
	//compute_h2b_computations();
	  start = std::clock();
	  compute_h2b_computations();
	  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	  mexPrintf("\tcompute_h2b_computations %f\n",duration);
  }
  if (penalty_weights[3]>0){
	//compute_h2c_computations();
	  start = std::clock();
	compute_h2c_computations();
	  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	  mexPrintf("\tcompute_h2c_computations %f\n",duration);
  }
  */
  energy=0.0;
  if (penalty_weights[0]>0){
	energy+=penalty_weights[0]*poisson_loglikelihood;
  }
  if (penalty_weights[1]>0){
	energy+=penalty_weights[1]*h1_uniform_spacing_value;
  }
  if (penalty_weights[2]>0){
	 energy+=penalty_weights[2]*h2b_pairwise_spacing_value;
  }
  if (penalty_weights[3]>0){
	energy+=penalty_weights[3]*h2c_pairwise_spacing_value;
  }


  return energy;
}

void SIMBA3D::simba3d::compute_gradient(){
	compute_poisson_loglikelihood_gradient();
	compute_h1_gradient();
	compute_h2b_gradient();
	compute_h2c_gradient();
}
void SIMBA3D::simba3d::compute_gradient(double *grad_x,double *grad_y,double *grad_z){
  //std::clock_t start;
  //double duration;
  if (penalty_weights[0]>0){
    //start = std::clock();
	compute_poisson_loglikelihood_gradient_0(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_poisson_loglikelihood_gradient %f\n",duration);
  }
  if (penalty_weights[1]>0){
    //start = std::clock();
	compute_h1_gradient(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h1_gradient %f\n",duration);
  }
  if (penalty_weights[2]>0){
    //start = std::clock();
	//compute_h2b_computations();
	compute_h2b_gradient(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h2b_gradient %f\n",duration);
  }
  if (penalty_weights[3]>0){
    //start = std::clock();
	//compute_h2c_computations();
	compute_h2c_gradient(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h2c_gradient %f\n",duration);
  }
}
void SIMBA3D::simba3d::compute_gradient(double number_of_threads,double *grad_x,double *grad_y,double *grad_z){
  //std::clock_t start;
  //double duration;
  if (penalty_weights[0]>0){
    //start = std::clock();
	compute_poisson_loglikelihood_gradient_0(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_poisson_loglikelihood_gradient %f\n",duration);
  }
  if (penalty_weights[1]>0){
    //start = std::clock();
	compute_h1_gradient(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h1_gradient %f\n",duration);
  }
  if (penalty_weights[2]>0){
    //start = std::clock();
	//compute_h2b_computations();
	compute_h2b_gradient(grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h2b_gradient %f\n",duration);
  }
  if (penalty_weights[3]>0){
    //start = std::clock();
	//compute_h2c_computations();
	compute_h2c_gradient(number_of_threads,grad_x,grad_y,grad_z);
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("\tcompute_h2c_gradient %f\n",duration);
  }
}
void SIMBA3D::simba3d::get_gradient(double *grad_x,double *grad_y,double *grad_z){
	for (int ii=0;ii<length;ii++){
		grad_x[ii]=x_gradient[ii];
		grad_y[ii]=y_gradient[ii];
		grad_z[ii]=z_gradient[ii];
	}
}
