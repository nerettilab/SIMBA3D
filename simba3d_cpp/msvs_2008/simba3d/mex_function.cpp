#include "mex.h"
#include <iostream>
#include <vector>
#include "simba3d.hpp"
#include <ctime>


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x,*y;
  size_t mrows,ncols;


  SIMBA3D::simba3d simba3d_instance=SIMBA3D::simba3d();      

  //mexPrintf("%d\n",nlhs);
  for (int ii=0;ii<nrhs;ii++){
	mrows = mxGetM(prhs[ii]);
	ncols = mxGetN(prhs[ii]);
	//mexPrintf("%d\t",mrows );
	//mexPrintf("%d\n",ncols );
  }
  int ii;
  // get penalty weights
  ii=0;
  mrows = mxGetM(prhs[ii]);
  ncols = mxGetN(prhs[ii]);
  simba3d_instance.set_penalty_weights(mrows*ncols,mxGetPr(prhs[ii]));
  // get penalty parameters
  ii=1;
  mrows = mxGetM(prhs[ii]);
  ncols = mxGetN(prhs[ii]);
  simba3d_instance.set_penalty_parameters(mrows*ncols,mxGetPr(prhs[ii]));  
  // get curve
  ii=2; 
  mrows = mxGetM(prhs[ii]);
  ncols = mxGetN(prhs[ii]);
  simba3d_instance.set_curve(mrows*ncols,mxGetPr(prhs[ii]),mxGetPr(prhs[ii+1]),mxGetPr(prhs[ii+2]));
  
  // get contact data
  ii=5; 
  mrows = mxGetM(prhs[ii]);
  ncols = mxGetN(prhs[ii]);
  simba3d_instance.set_contact(mrows*ncols,mxGetPr(prhs[ii]),mxGetPr(prhs[ii+1]),mxGetPr(prhs[ii+2]));
  
  std::clock_t start;
  double duration;
  double * energy;

  start = std::clock();
  // perform initial computations
  if (simba3d_instance.get_length()>1000){
	simba3d_instance.compute_pairwise_computations_0(4);
  } else {// otherwise it is not worth the overhead for doing the threading
	simba3d_instance.compute_pairwise_computations();
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  //mexPrintf("compute_pairwise_computations_0 %f\n",duration);

  start = std::clock();
  // perform initial computations
  simba3d_instance.compute_adjacent_computations();
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  //mexPrintf("compute_adjacent_computations %f\n",duration);

  plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  energy=mxGetPr(plhs[0]);

  if (nlhs==1){
	  simba3d_instance.compute_poisson_loglikelihood();
  } else {
  }

      

  simba3d_instance.compute_h2b_computations();
  simba3d_instance.compute_h2c_computations();

  //simba3d_instance.create_index_mutex();
  if (nlhs==4){
    //simba3d_instance.create_index_mutex();
	plhs[1] = mxCreateDoubleMatrix((mwSize)simba3d_instance.get_length(), (mwSize)1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)simba3d_instance.get_length(), (mwSize)1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((mwSize)simba3d_instance.get_length(), (mwSize)1, mxREAL);
	double *grad_x=mxGetPr(plhs[1]);
	double *grad_y=mxGetPr(plhs[2]);
	double *grad_z=mxGetPr(plhs[3]);
	memset(grad_x, 0, sizeof(grad_x));
	memset(grad_y, 0, sizeof(grad_y));
	memset(grad_z, 0, sizeof(grad_z));
	
	start = std::clock();
	simba3d_instance.compute_gradient(4,grad_x,grad_y,grad_z);
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//mexPrintf("compute_gradient %f\n",duration);
	
	/*
	simba3d_instance.clear_gradient();
	simba3d_instance.compute_gradient();
	simba3d_instance.get_gradient(grad_x,grad_y,grad_z);
	*/
  }
   simba3d_instance.compute_energy();
  *energy  = simba3d_instance.get_engery();
  /* Check for proper number of arguments.
  if(nrhs!=1) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:invalidNumInputs",
            "One input required.");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:maxlhs",
            "Too many output arguments.");
  }
   */
  /* The input must be a noncomplex scalar double.
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(toVector(prhs[0])) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Input must be a noncomplex scalar double.");
  }
  */
  //simba3d_instance.compute_pairwise_computations();
  //simba3d_instance.compute_adjacent_computations();
  //simba3d_instance.clear_gradient();
  //simba3d_instance.print_gradient();
  //simba3d_instance.compute_poisson_loglikelihood();
  //simba3d_instance.compute_poisson_loglikelihood_gradient();
  //simba3d_instance.compute_h1_gradient();
  //simba3d_instance.print_gradient();
  //simba3d_instance.compute_h2b_computations();
  //simba3d_instance.compute_h2b_gradient();
  //simba3d_instance.print_gradient();
  //simba3d_instance.compute_h2c_computations();
  //simba3d_instance.compute_h2c_gradient();
  //simba3d_instance.print_gradient();
  /* Create matrix for the return argument. */
  //plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
  
  /* Assign pointers to each input and output. */
  //x = mxGetPr(prhs[0]);
  //y = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */

}
