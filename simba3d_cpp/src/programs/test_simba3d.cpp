/**
* @file test_simba3d.cpp
* @brief Try and feed tena something, but she will not eat
*
* 
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "pthread.h"

#include "simba3d.hpp"

int main(int argc, char *argv[]){

  // create instance class object
  SIMBA3D::simba3d simba3d_instance=SIMBA3D::simba3d();      

  // input the penalty weights  
  std::vector<double> penalty_weights;
  // data term weight
  penalty_weights.push_back(1.0);
  // h1 uniform spacing penalty
  penalty_weights.push_back(1.0);
  // h2 smoothing penalty
  penalty_weights.push_back(0.0);  
  // h3 population prior penalty
  penalty_weights.push_back(1.0);  
  // h4 lamina penalty
  penalty_weights.push_back(0.0);    
  // h5 repulsion penalty
  penalty_weights.push_back(0.0);    
  // set the term weights
  simba3d_instance.set_penalty_weights(penalty_weights);

  // make an example curve
  std::vector<double> x,y,z;
  //
  x.push_back(0.0);
  x.push_back(2.0);
  x.push_back(3.0);
  x.push_back(0.1);
  //
  y.push_back(0.0);
  y.push_back(1.0);
  y.push_back(2.0);
  y.push_back(0.1);
  //
  z.push_back(0.0);
  z.push_back(1.0);
  z.push_back(2.0);
  z.push_back(0.1);
  /*
  for (int ii=0;ii<10000;ii++){
	x.push_back(2.0);
	y.push_back(3.0);
	z.push_back(0.0);
  }
  */
  pthread_mutex_t index_locks= PTHREAD_MUTEX_INITIALIZER;
  
  std::cout << pthread_mutex_lock(&index_locks) <<"\n"; 
  std::cout << pthread_mutex_unlock(&index_locks) <<"\n"; 
 
  simba3d_instance.set_curve(x,y,z);

  std::vector<double> row,col,value;
  row.push_back(0);
  col.push_back(1);
  value.push_back(15.0);
  /*
  row.push_back(0);
  col.push_back(1);
  value.push_back(105.0);
  row.push_back(0);
  col.push_back(2);
  value.push_back(0.0);
  row.push_back(0);
  col.push_back(3);
  value.push_back(0.0);
  row.push_back(1);
  col.push_back(2);
  value.push_back(0.0);
  row.push_back(1);
  col.push_back(3);
  value.push_back(0.0);
  row.push_back(2);
  col.push_back(3);
  value.push_back(0.0);
  */
  simba3d_instance.set_contact(row,col,value);
  //
  simba3d_instance.compute_pairwise_computations();

  simba3d_instance.compute_adjacent_computations();
  simba3d_instance.compute_h2b_computations();
  simba3d_instance.compute_h2c_computations();
  simba3d_instance.clear_gradient();
  //simba3d_instance.print_gradient();
  simba3d_instance.compute_gradient();
  double e=simba3d_instance.compute_energy();
  
  std::cout << e << "\n";
  simba3d_instance.print_gradient();

  return 0;
}
