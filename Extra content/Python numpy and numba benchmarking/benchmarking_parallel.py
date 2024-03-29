# -*- coding: utf-8 -*-
"""
Benchmarking the performance of the various versions of Solve_HH_Problem 
used to solve the stochastic growth model.

This file runs the loop over states in parallel. Due to the model features, 
only the base for-loop (version 1) function runs faster in parallel.

Using 1000 iterations takes approximately 20 minutes to run.
"""

import numpy as np
from datetime import datetime

# Import parameter values and struct.
from main import (alpha, beta, delta, k_values, number_of_k_values, z_values, z_probs, 
                  number_of_z_values, Parameters)
params = Parameters(alpha, beta, delta, k_values, z_values, z_probs)

# Assign additional parameter values.
number_of_iterations = 1000

# Initialize Value Function and Policy Function (as arrays).
Value_Function = np.zeros((number_of_iterations, number_of_k_values, number_of_z_values))
Policy_Function = np.zeros((number_of_iterations, number_of_k_values, number_of_z_values))

# Load the various versions of the key value function iteration function.
from custom_functions import (Solve_HH_Problem_v1, Solve_HH_Problem_v2,
                              Solve_HH_Problem_v3, Solve_HH_Problem_v4)

functions = [("base", Solve_HH_Problem_v1),
             ("numpy", Solve_HH_Problem_v2),
             ("numba precompiled", Solve_HH_Problem_v3),
             ("numba precompiled + numpy", Solve_HH_Problem_v4)]


# Solve model with each version of the function.
for (version, Solve_HH_Problem) in functions:
    
    start = datetime.now()
    
    import os
    from joblib import Parallel, delayed
    num_threads = os.cpu_count()

    if __name__ == "__main__":
        
        for iteration in range(1, number_of_iterations):
            
            # Get list of Value Function and Policy Function values for each starting state (kt0, zt).
            processed_list = Parallel(n_jobs=num_threads)(delayed(Solve_HH_Problem)
                                                        (Value_Function[iteration-1, : :], kt0, zt, params)
                                                        for kt0 in enumerate(k_values) for zt in enumerate(z_values))
            
            # Reshape to arrays in the same shape as the Value Function and Policy Function.
            new_values = np.reshape(processed_list, (number_of_k_values, number_of_z_values, 2))
            
            # Update values.
            Value_Function[iteration,:,:] = new_values[:, :, 0]
            Policy_Function[iteration,:,:] = new_values[:, :, 1]
    
    end = datetime.now()
    
    # Display results.
    print(version)
    print("   start: ", start.strftime("%H:%M:%S"))
    print("   stop: ", end.strftime("%H:%M:%S"))
    print("   seconds: ", (end-start).seconds)
    print("   minutes: ", round((end-start).seconds/60,2))
    print("\n")