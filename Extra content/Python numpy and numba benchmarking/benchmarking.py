# -*- coding: utf-8 -*-
"""
Benchmarking the performance of the various versions of Iterate_Value_Function 
used to solve the stochastic growth model.

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
    
    # Part 1 - implement value function iteration without parallel processing.
    start = datetime.now()

    # Perform value function iteration.
    for iteration in range(1,number_of_iterations):
        
        # Loop over all possible starting states.
        for kt0 in enumerate(k_values):
            for zt in enumerate(z_values):
                
                # Solve Value Function and Policy Function and update values
                V, g = Solve_HH_Problem(Value_Function[iteration-1, :, :], kt0, zt, params)
                
                Value_Function[iteration, kt0[0], zt[0]]  = V
                Policy_Function[iteration, kt0[0], zt[0]] = g
                
    end = datetime.now()
    
    # Display results.
    print(version)
    print("Single-core implementation:")
    print("   start: ", start.strftime("%H:%M:%S"))
    print("   stop: ", end.strftime("%H:%M:%S"))
    print("   seconds: ", (end-start).seconds)
    print("   minutes: ", (round((end-start).seconds/60,2)))