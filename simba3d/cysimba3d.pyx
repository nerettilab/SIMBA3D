# distutils: language = c++
# distutils: sources = ./simba3d_cpp/src/simba3d.cpp
import numpy as np
import time
from cpython cimport array
cimport numpy as np

cdef extern from "../simba3d_cpp/src/simba3d.hpp" namespace "SIMBA3D":
    cdef cppclass simba3d:
        simba3d() except +
        void set_poisson_weight(double value) except +
        void set_h1_weight(double value) except +
        void set_h2a_weight(double value) except +
        void set_h2b_weight(double value) except +
        void set_h2c_weight(double value) except +
        void set_number_of_subtasks(int number)  except +
        void set_contact(int input_number_of_contacts,double * row_index, double * column_index,double * value ) except +
        int get_length() except +
        void set_curve(int input_length,double * input_x,double * input_y,double * input_z) except +
        void initialize_pairwise_computations() except +
        void compute_pairwise_computations() except +
        void compute_pairwise_computations(int number_of_threads) except +
        void compute_pairwise_computations_0(int number_of_threads) except +
        void compute_adjacent_computations() except +
        void compute_poisson_loglikelihood() except +
        void compute_h2b_computations() except +
        void compute_h2c_computations() except +
        double compute_energy() except +
        double get_engery() except +
        void compute_gradient(double number_of_threads,double *grad_x,double *grad_y,double *grad_z) except +
        double compute_poisson_loglikelihood_gradient_0(double *grad_x,double *grad_y,double *grad_z) except +
        double compute_poisson_loglikelihood_gradient_1(double *grad_x,double *grad_y,double *grad_z) except +
        double compute_h1_gradient(double *grad_x,double *grad_y,double *grad_z) except +
        double compute_h2a_gradient(double *grad_x,double *grad_y,double *grad_z) except +
        double compute_h2b_gradient(double *grad_x,double *grad_y,double *grad_z) except +
        double compute_h2c_gradient(int number_of_threads, double *grad_x,double *grad_y,double *grad_z) except +
        void print_gradient() except +

cdef class gradient_manager:
    cdef simba3d c_obj  # hold a C++ instance which we're wrapping
    cdef int number_of_threads
    cdef list functions
    cdef list shared_computations   
    def __cinit__(self,term_parameters=None,number_of_threads=None):
        self.c_obj =  simba3d()
        if number_of_threads is None:
            self.number_of_threads=1
        else:
            self.number_of_threads=number_of_threads
        print(self.number_of_threads)
        self.functions=[]
        self.shared_computations=[]
        if term_parameters is None:
            term_parameters={}
        print('Reading Term Weights')
        pairwise_computations_needed=False
        adjacent_computations_needed=False
        h2b_computations_needed=False
        h2c_computations_needed=False
        print(term_parameters.keys())
        if "h1_weight" in term_parameters:
            self.c_obj.set_h1_weight(term_parameters["h1_weight"])   
            if term_parameters["h1_weight"]>0:
                #pairwise_computations_needed=True 
                adjacent_computations_needed=True 
                self.functions.append(self.compute_h1_gradient)    
        if "h2a_weight" in term_parameters:
            self.c_obj.set_h2a_weight(term_parameters["h2a_weight"]) 
            if term_parameters["h2a_weight"]>0:
                #pairwise_computations_needed=True 
                adjacent_computations_needed=True
                self.functions.append(self.compute_h2a_gradient)  
        if "h2b_weight" in term_parameters:
            self.c_obj.set_h2b_weight(term_parameters["h2b_weight"])   
            if term_parameters["h2b_weight"]>0:
                pairwise_computations_needed=True 
                h2b_computations_needed=True
                self.functions.append(self.compute_h2b_gradient)   
        if "h2c_weight" in term_parameters:
            self.c_obj.set_h2c_weight(term_parameters["h2c_weight"]) 
            if term_parameters["h2c_weight"]>0:
                pairwise_computations_needed=True   
                h2c_computations_needed=True
                self.functions.append(self.compute_h2c_gradient)   
        if "poisson_weight" in term_parameters:
            self.c_obj.set_poisson_weight(term_parameters["poisson_weight"])
            if term_parameters["poisson_weight"]>0:
                if pairwise_computations_needed:
                    self.functions.append(self.compute_poisson_loglikelihood_gradient_0)
                else:
                    self.functions.append(self.compute_poisson_loglikelihood_gradient_1)
        self.shared_computations.append( self.initialize_pairwise_computations)
        #self.shared_computations.append()
        if pairwise_computations_needed:
            self.shared_computations.append(self.run_pairwise_computations)
        if adjacent_computations_needed:
            self.shared_computations.append(self.run_adjacent_computations)
        if h2b_computations_needed:
            self.shared_computations.append(self.run_h2b_computations)
        if h2c_computations_needed:
            self.shared_computations.append(self.run_h2b_computations)
            self.shared_computations.append(self.run_h2c_computations)
    def initialize_pairwise_computations(self):
        return self.c_obj.initialize_pairwise_computations()
    def compute_poisson_loglikelihood_gradient_0(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        return self.c_obj.compute_poisson_loglikelihood_gradient_0(&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def compute_poisson_loglikelihood_gradient_1(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        return self.c_obj.compute_poisson_loglikelihood_gradient_1(&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def compute_h1_gradient(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        return self.c_obj.compute_h1_gradient(&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def compute_h2a_gradient(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        return self.c_obj.compute_h2a_gradient(&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def compute_h2b_gradient(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        return self.c_obj.compute_h2b_gradient(&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def compute_h2c_gradient(self,                        
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None):
        #return self.c_obj.compute_h2c_gradient(1,&x_gradient[0],&y_gradient[0],&z_gradient[0])
        return self.c_obj.compute_h2c_gradient(self.number_of_threads,&x_gradient[0],&y_gradient[0],&z_gradient[0])
    def run_pairwise_computations(self):
        '''
        '''
        self.c_obj.compute_pairwise_computations_0(self.number_of_threads);
        '''
        if  self.c_obj.get_length()>1000:
            self.c_obj.compute_pairwise_computations_0(self.number_of_threads);
        else:# otherwise it is not worth the overhead for doing the threading
            self.c_obj.compute_pairwise_computations();        
        '''
    def run_adjacent_computations(self):
        '''
        '''
        self.c_obj.compute_adjacent_computations()
    def run_h2b_computations(self):
        '''
        '''
        self.c_obj.compute_h2b_computations()
    def run_h2c_computations(self):
        '''
        '''
        self.c_obj.compute_h2c_computations()
    def set_contact_data(self,contact_row_index,contact_col_index,contact_value):
        cdef array.array r=array.array('d',contact_row_index)
        cdef array.array c=array.array('d',contact_col_index)
        cdef array.array v=array.array('d',contact_value)
        self.c_obj.set_contact(len(contact_row_index),r.data.as_doubles,c.data.as_doubles,v.data.as_doubles)
    def compute_gradient(self,
                        np.ndarray[double,ndim=1,mode="c"] x not None,
                        np.ndarray[double,ndim=1,mode="c"] y not None,
                        np.ndarray[double,ndim=1,mode="c"] z not None,
                        np.ndarray[double,ndim=1,mode="c"] x_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] y_gradient not None,
                        np.ndarray[double,ndim=1,mode="c"] z_gradient not None
                        ):
        self.c_obj.set_curve(len(x),&x[0],&y[0],&z[0])
        # Loop through the list of shared calculations
        self.c_obj.set_number_of_subtasks(self.number_of_threads)
        for shared_computation in self.shared_computations:
            #print(shared_computation)
            shared_computation()

        #self.c_obj.compute_pairwise_computations_0(4)
        #self.c_obj.compute_adjacent_computations();
        #self.c_obj.compute_poisson_loglikelihood();
        #self.c_obj.compute_h2b_computations();
        #self.c_obj.compute_h2c_computations();        
        # loop through the non-trivial terms
        e=0.0
        for functional in self.functions:
            #print(functional)
            e+=functional(x_gradient,y_gradient,z_gradient)
        #self.c_obj.compute_gradient(self.number_of_threads,&x_gradient[0],&y_gradient[0],&z_gradient[0])
       # self.c_obj.compute_energy()
        return e
    def get_engery(self):
        return self.c_obj.get_engery()
   