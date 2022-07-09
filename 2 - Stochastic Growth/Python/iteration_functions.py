# -*- coding: utf-8 -*-

"""
Function to perform an iteration of value function iteration.
Four different versions with increasing levels of performance.
Each version has identical inputs and outputs.

Input:
    - Current value function
    - Starting capital and productivity level 
    
Output:
    - Next value function value
    - Optimal choice of next period capital
    
Version
    v1: using only base functions + for-loops
    v2: using Numpy arrays / broadcast instead of for-loop
    v3: precompiled with numba and using for-loops
    v4: precompiled with numba and using numpy arrays instead of for-loop
"""

import numpy as np
from math import log

# -----------------------------------------------------------------------------------------------------
# Version 1 - using only base functions + for-loops.
# -----------------------------------------------------------------------------------------------------
def Iterate_Value_Function_v1(Previous_Value_Function, kt0, zt, params):
    
    # Unpack utility parameters and grids.
    alpha,beta,delta = params.alpha, params.beta, params.delta
    k_values = params.k_values
    z_probs = params.z_probs
    
    # Variables to store candidate optimal values for Value Function and Policy Function.
    v_max = -float('inf')
    kt1_optimal = 0.0
        
    # Search over possible next period capital choices.
    for kt1 in enumerate(k_values):
        
        # Calculate value function for given choice of next period capital.
        new_v_max = ( log( (zt[1]*(kt0[1] ** alpha)) + (1-delta)*kt0[1] - kt1[1]) 
                     + beta*np.dot(Previous_Value_Function[kt1[0],:], z_probs[zt[0],:]))
        
        # Check if this capital choice gives highest Value Function value.
        if new_v_max > v_max:
        
            # Update candidate values.
            v_max = new_v_max
            kt1_optimal = kt1[1]
    
    return v_max, kt1_optimal


# -----------------------------------------------------------------------------------------------------
# Version 2 - using Numpy arrays / broadcast instead of for-loop.
# -----------------------------------------------------------------------------------------------------
def Iterate_Value_Function_v2(Previous_Value_Function, kt0, zt, params):
    
    # Unpack utility parameters and grids.
    alpha,beta,delta = params.alpha, params.beta, params.delta
    k_values = params.k_values
    z_probs = params.z_probs
    
    # Calculate array of value function values for all k_values.
    V_max_values = (np.log( zt[1]*(kt0[1]**alpha) + (1-delta)*kt0[1] - k_values) 
                    + beta*np.dot(Previous_Value_Function, z_probs[zt[0],:]))
    
    # Get the index for the optimal capital choice.
    kt1_index_optimal = np.argmax(V_max_values)
    
    # Get Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_index_optimal]
    v_max = V_max_values[kt1_index_optimal]
    
    return v_max, kt1_optimal


# -----------------------------------------------------------------------------------------------------
# Version 3 - precompiled with numba and using for-loops.
# -----------------------------------------------------------------------------------------------------
from numba import jit

@jit(nopython=True)    
def Iterate_Value_Function_v3(Previous_Value_Function, kt0, zt, params):
    
    # Unpack utility parameters and grids.
    alpha,beta,delta = params.alpha, params.beta, params.delta
    k_values = params.k_values
    z_probs = params.z_probs
    
    # Variables to store candidate optimal values for Value Function and Policy Function.
    v_max = 0.0
    new_v_max = 0.0
    kt1_optimal = 0.0
        
    # Search over possible next period capital choices
    for kt1 in enumerate(k_values):
        
        # Calculate value function for choice of next period capital.
        new_v_max = ( log((zt[1]*(kt0[1] ** alpha) + (1-delta)*kt0[1] - kt1[1])) 
                     + beta*np.dot(Previous_Value_Function[kt1[0],:], z_probs[zt[0],:]))
        
        # Check if this capital choice gives highest Value Function value.
        if new_v_max > v_max:
        
            # Update candidate values.
            v_max = new_v_max
            kt1_optimal = kt1[1]
    
    return v_max, kt1_optimal 


# -----------------------------------------------------------------------------------------------------
# Version 4 - precompiled with numba and using numpy arrays instead of for-loop.
# -----------------------------------------------------------------------------------------------------    
from numba import jit

@jit(nopython=True)                      
def Iterate_Value_Function_v4(Previous_Value_Function, kt0, zt, params):
    
    # Unpack utility parameters and grids.
    alpha,beta,delta = params.alpha, params.beta, params.delta
    k_values = params.k_values
    z_probs = params.z_probs
    
    # Calculate array of value function values for all k_values.
    V_max_values = (np.log( zt[1]*(kt0[1]**alpha) + (1-delta)*kt0[1] - k_values) 
                    + beta*np.dot(Previous_Value_Function, z_probs[zt[0],:]))
    
    # Get the index for the optimal capital choice.
    kt1_ind_optimal = np.argmax(V_max_values)
    
    # Get Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_ind_optimal]
    v_max = V_max_values[kt1_ind_optimal]
    
    return v_max, kt1_optimal