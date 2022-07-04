# -*- coding: utf-8 -*-

"""
Function to perform an iteration of value function iteration.
Four different versions with increasing levels of performance.
Each version has idential inputs and outputs.

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

# Import parameter values.
from Python_implementation import (alpha, beta, delta, k_values, z_probs)
import numpy as np
from math import log

"""
Version 1 - using only base functions + for-loops.
"""

def Iterate_Value_Function_v1(Previous_Value_Function, kt0, zt):
    
    # Variables to store candidate optimal values for Value Function and Policy Function.
    v_max = -float('inf')
    kt1_optimal = 0.0
        
    # Search over possible next period capital choices
    for kt1 in enumerate(k_values):
        
        # Calculate value function for choice of next period capital.
        new_v_max = ( log( (zt[1]*(kt0[1] ** alpha)) + (1-delta)*kt0[1] - kt1[1]) 
                     + beta*np.dot(Previous_Value_Function[kt1[0],:], z_probs[zt[0],:]))
        
        # Check if this capital choice gives highest Value Function value.
        if new_v_max > v_max:
        
            # Update candidate values.
            v_max = new_v_max
            kt1_optimal = kt1[1]
    
    return v_max, kt1_optimal


"""
Version 2 - using Numpy arrays / broadcast instead of for-loop.
"""

def Iterate_Value_Function_v2(Previous_Value_Function, kt0, zt):
    
    # Calculate array of value function values for all k_values.
    V_max_values = (np.log( zt[1]*(kt0[1]**alpha) + (1-delta)*kt0[1] - k_values) 
                    + beta*np.dot(Previous_Value_Function, z_probs[zt[0],:]))
    
    # Get the index for the optimal capital choice.
    kt1_index_optimal = np.argmax(V_max_values)
    
    # Get Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_index_optimal]
    v_max = V_max_values[kt1_index_optimal]
    
    return v_max, kt1_optimal


"""
Version 3 - precompiled with numba and using for-loops.
"""           
      
from numba import jit

@jit(nopython=True)    
def Iterate_Value_Function_v3(Previous_Value_Function, kt0, zt):
    
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


"""
Version 4 - precompiled with numba and using numpy arrays instead of for-loop.
"""
            
from numba import jit

@jit(nopython=True)                      
def Iterate_Value_Function_v4(Previous_Value_Function, kt0, zt):
    
    # Calculate array of value function values for all k_values.
    V_max_values = (np.log( zt[1]*(kt0[1]**alpha) + (1-delta)*kt0[1] - k_values) 
                    + beta*np.dot(Previous_Value_Function, z_probs[zt[0],:]))
    
    # Get the index for the optimal capital choice.
    kt1_ind_optimal = np.argmax(V_max_values)
    
    # Get Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_ind_optimal]
    v_max = V_max_values[kt1_ind_optimal]
    
    return v_max, kt1_optimal

  
"""
Alternative to importing parameters: 
    re-define parameters here. May be needed to use parallel implementation.
"""
"""
# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 0.012

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
"""