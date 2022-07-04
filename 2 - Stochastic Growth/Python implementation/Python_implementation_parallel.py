# -*- coding: utf-8 -*-
"""
Benchmarking the performance of the various versions of Iterate_Value_Function 
used to solve the stochastic growth model.

This file runs the loop over states in parallel. Due to the model features, 
only the base for-loop (version 1) funtion runs faster in parallel.

Using 1000 iterations takes approximately 20 minutes to run.
"""

import numpy as np
from datetime import datetime

# Parallel implementation requires definition of these variables in this file.
# It may also require definition ofthese variables in the iteration_function module.

# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 0.012
number_of_iterations = 10

# Calculate the steady-state level of capital.
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ** (1/(alpha-1))  

# Create a range of capital values around steady-state (+/- 2%).
number_of_k_values = 201
k_low_pct = 0.98
k_high_pct = 1.02
k_values = np.linspace(k_low_pct*k_steady, k_high_pct*k_steady, num=number_of_k_values)

# Get productivity levels and transition probabilities.
z_probs = np.genfromtxt("Inputs\z_probs.csv", delimiter=",")
z_values = np.genfromtxt("Inputs\z_values.csv", delimiter=",")
number_of_z_values = len(z_values)

# Initialize Value Function and Policy Function (as arrays).
Value_Function = np.zeros((number_of_iterations, number_of_k_values, number_of_z_values))
Policy_Function = np.zeros((number_of_iterations, number_of_k_values, number_of_z_values))


# Load the various versions of the key value function iteration function.
import iteration_functions
from iteration_functions import (Iterate_Value_Function_v1, Iterate_Value_Function_v2,
                                 Iterate_Value_Function_v3, Iterate_Value_Function_v4)

functions = [("base", Iterate_Value_Function_v1),
             ("numpy", Iterate_Value_Function_v2),
             ("numba precompiled", Iterate_Value_Function_v3),
             ("numba precompiled + numpy", Iterate_Value_Function_v4)]


# Solve model with each version of the function.
for (version, Iterate_Value_Function) in functions:
    
    start = datetime.now()
    
    import os
    from joblib import Parallel, delayed
    num_threads = os.cpu_count()

    if __name__ == "__main__":
        
        for iteration in range(1, number_of_iterations):
            
            # Get list of Value Function and Policy Function values for each starting state (kt0, zt).
            processed_list = Parallel(n_jobs=num_threads)(delayed(Iterate_Value_Function)
                                                        (Value_Function[iteration-1,::], kt0, zt)
                                                        for kt0 in enumerate(k_values) for zt in enumerate(z_values))
            
            # Reshape to arrays in the same shape as the Value Function and Policy Function.
            new_values = np.reshape(processed_list, (number_of_k_values, number_of_z_values, 2))
            
            # Update values.
            Value_Function[iteration,:,:] = new_values[:,:,0]
            Policy_Function[iteration,:,:] = new_values[:,:,1]
    
    end = datetime.now()
    
    # Display results.
    print(version)
    print("   start: ", start.strftime("%H:%M:%S"))
    print("   stop: ", end.strftime("%H:%M:%S"))
    print("   seconds: ", (end-start).seconds)
    print("   minutes: ", round((end-start).seconds/60,2))
    print("\n")