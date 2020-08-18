/**
* @file simba3d.hpp
* @brief a brief description for  the template name header file 
*
* If you want, you can say more about the contents of the header. 
*
*
*/

#ifndef SIMBA3D_H_
#define SIMBA3D_H_

#include "pthread.h"
#include <string>
#include <vector>
#include <cmath>
#include <iostream>



namespace SIMBA3D{
            
void split_pairwise_distance_computation(
  double* pairwise_distances,
  int number_of_pairwise_distances,
  int division,
  int length
  );

void triu_ij_to_ind(int i , int j,int & ind);
void triu_ind_to_ij(int &i , int & j,int ind);

void run_pairwise_computations( 
                int length,
                double *x,
                double *y,
                double *z,
                int number_of_pairs,
                double *pairwise_difference_x,
                double *pairwise_difference_y,
                double *pairwise_difference_z,
                double *pairwise_distance,
                double *series
                );

void compute_pairwise_difference(
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
);
  /**
  \brief version number
  */
  const char VERSION[]="00.00.00";

  /**
  \brief store shared data for computing simba3d gradients
  

  [example_feed_tena_dinner](../../src/programs/test_simba3d.cpp)
  */
  class simba3d {
    private:
      // inputs......................................................
      std::vector<double> penalty_weights;
      int penalty_inventory_length;     
  	  std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      int number_of_contact_pairs;
      std::vector<int> contact_row_index;
      std::vector<int> contact_column_index;
      std::vector<double> contact_value;
      double contact_sum;
      /**< \brief upper triangular differences*/   
      int length;
      /**< \brief the number of  node points for the 3D structure*/  

      // pairwise calculations .......................................
      int number_of_pairs;
      double sum_pairwise_distance;
      double average_pairwise_distance;
      std::vector<double> x_differences;
      /**< \brief upper triangular differences*/
      std::vector<double> y_differences;
      /**< \brief upper triangular differences*/
      std::vector<double> z_differences;
      /**< \brief upper triangular differences*/
      std::vector<double> distances;
      /**< \brief upper triangular distances*/
      // data term
      double poisson_a_parameter;
      double poisson_b_parameter;
      double poisson_loglikelihood;
      // adjacent calculations  ......................................
      double sum_adjacent_distance;
      double average_adjacent_distance;
      double sum_of_squared_adjacent_distance;
      // h1 scale invariant uniform spacing penalty
      double h1_uniform_spacing_value;
      double h2a_value;
      // h2 (b) scale invariant pairwise spacing penalty
      double h2b_pairwise_spacing_value;
      double alpha_h2b;
      std::vector<double> F_h2b;
      std::vector<double> H_h2b;
      std::vector<double> G_h2b;
      // h2 (c) scale invariant pairwise spacing penalty with threshold
      double h2c_pairwise_spacing_value;
      double h2c_radius;
      std::vector<int> h2c_index_within_threshold;
      // multi threading
      int number_of_subtasks;
      std::vector<double> split_sum_pairwise_distance;
	    std::vector<pthread_mutex_t> index_locks;
     // outputs .....................................................
  		std::vector<double> x_gradient;
  		std::vector<double> y_gradient;
  		std::vector<double> z_gradient;      
    		double energy;

        /**< \brief What tena has for dinner.*/
    public:

    /**
    \brief Initialize an instance without arguments (using default dinner).
    */  
	simba3d();
    /**
    \brief Destruct instance of simba3d.

    */  
	~simba3d();
    /**
    \brief set the 3D curve
    */
  void initialize_penalty_weights();       
	void set_penalty_weights(
	   const std::vector<double> & input_penalty_weights
	  /**< \brief a vectorized 3d curve.
		

	  */
	);
	void set_penalty_weights(
	   int input_number_of_penalties,
	   double * input_penalty_weights
	  /**< \brief a vectorized 3d curve.
		

	  */
	); 
  void set_poisson_weight(double value){penalty_weights[0]=value;};   
  void set_h1_weight(double value){penalty_weights[1]=value;};   
  void set_h2a_weight(double value){penalty_weights[4]=value;};
  void set_h2b_weight(double value){penalty_weights[2]=value;};
  void set_h2c_weight(double value){penalty_weights[3]=value;};
	void set_penalty_parameters(
	   int input_number_of_penalties,
	   double * input_penalty_parameters
	  /**< \brief a vectorized 3d curve.
		

	  */
	); 
	void set_curve(
    const std::vector<double> & input_x,
    /**< \brief x coordinate*/
    const std::vector<double> & input_y,
    /**< \brief y coordinate*/
    const std::vector<double> & input_z
    /**< \brief z coordinate*/
   );
	void set_curve(
	int input_length,
    double * input_x,
    /**< \brief x coordinate*/
    double * input_y,
    /**< \brief y coordinate*/
    double * input_z
    /**< \brief z coordinate*/
   );	
  void set_pairwise(
    int input_length,
    /**< \brief the number of nodes 
      note that the size of the arrays passed in should be length*(length-1)/2
    */
    double * input_x_differences,
    /**< \brief upper triangular differences*/
    double * input_y_differences,
    /**< \brief upper triangular differences*/
    double * input_z_differences,
   /**< \brief upper triangular differences*/   
    double * input_distances
   /**< \brief upper triangular distances*/
  );
   void set_contact(
    const std::vector<double> & row_index,
    /**< \brief x coordinate*/
    const std::vector<double> & column_index,
    /**< \brief y coordinate*/
    const std::vector<double> & value
    /**< \brief z coordinate*/
  );   
   void set_contact(
    int input_number_of_contacts,
    double * row_index,
    /**< \brief x coordinate*/
    double * column_index,
    /**< \brief y coordinate*/
    double * value
    /**< \brief z coordinate*/
  );   
  void initialize_pairwise_computations();
  void compute_pairwise_computations();
  void compute_pairwise_computations(int number_of_threads);
  void compute_pairwise_computations_0(int number_of_threads);
  void set_number_of_subtasks(int number);
  void split_compute_pairwise_computations(int subtask_id);
  void combine_split_pairwise_sums();
  void compute_adjacent_computations();
  void compute_adjacent_computations_recompute();
  void clear_gradient();
  void clear_gradient(double *grad_x,double *grad_y,double *grad_z);
  void print_gradient();
  void compute_poisson_loglikelihood();
  void compute_poisson_loglikelihood_gradient();
  double compute_poisson_loglikelihood_gradient_0(double *grad_x,double *grad_y,double *grad_z);
  double compute_poisson_loglikelihood_gradient_1(double *grad_x,double *grad_y,double *grad_z);
  void compute_poisson_loglikelihood_gradient(double *grad_x,double *grad_y,double *grad_z);
  double compute_poisson_loglikelihood_gradient(int number_of_threads,double *grad_x,double *grad_y,double *grad_z);
  void split_compute_poisson_loglikelihood_gradient(int subtask_id,double *grad_x,double *grad_y,double *grad_z);
  void create_index_mutex();
  void compute_h1_gradient();
  double compute_h1_gradient(double *grad_x,double *grad_y,double *grad_z);
  double compute_h2a_gradient();
  double compute_h2a_gradient(double *grad_x,double *grad_y,double *grad_z);
  void compute_h2b_computations();
  void compute_h2b_gradient();
  double compute_h2b_gradient(double *grad_x,double *grad_y,double *grad_z);
  void compute_h2c_computations();
  void compute_h2c_gradient();
  void compute_h2c_gradient(double *grad_x,double *grad_y,double *grad_z);
  double compute_h2c_gradient(int number_of_threads,double *grad_x,double *grad_y,double *grad_z);
  void split_compute_h2c_gradient(int subtask_id,double *grad_x,double *grad_y,double *grad_z);
  int get_number_of_contact_pairs(){return  number_of_contact_pairs;}
  int get_contact_pairs(unsigned ind,int & r, int &c, int &v);
  int get_length(){return  length;}
  double compute_energy();
  double get_engery(){return energy;}
  void compute_gradient();
  void compute_gradient(double *grad_x,double *grad_y,double *grad_z);
  void compute_gradient(double number_of_threads,double *grad_x,double *grad_y,double *grad_z);
  void get_gradient(double *grad_x,double *grad_y,double *grad_z);
};



}
typedef struct {  int subtask_id;
          SIMBA3D::simba3d * instance;
		  double *grad_x;
		  double *grad_y;
		  double *grad_z;
	} thread_struct;

void *split_compute_pairwise_computations_helper(
	int subtask_id,
	SIMBA3D::simba3d * instance
	);     
void *split_compute_poisson_loglikelihood_gradient_helper(
	int subtask_id,
	SIMBA3D::simba3d * instance
	);    
void *split_compute_h2c_gradient_helper(
	int subtask_id,
	SIMBA3D::simba3d * instance
	);    


#endif  // SIMBA3D_H_
