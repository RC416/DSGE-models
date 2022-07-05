# -*- coding: utf-8 -*-
"""
Stochastic Growth model implemented in Python.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Create function to solve Value Function for given starting states.
    3 - Perform value function iteration.
    4 - Plot results.
"""

import numpy as np

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------
# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 0.012
number_of_iterations = 1000

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


# Store utility parameters and capital/productivity grids in a struct for passing to a function.
from typing import NamedTuple # could also use dataclass from dataclasses module

class Parameters(NamedTuple):
    alpha: float
    beta:  float
    delta: float
    k_values: float
    z_values: float
    z_probs:  float
    
params = Parameters(alpha, beta, delta, k_values, z_values, z_probs)

# -----------------------------------------------------------------------------------------------------
# 2 - Create function to solve Value Function for given starting states.
# -----------------------------------------------------------------------------------------------------
def Iterate_Value_Function(Previous_Value_Function, kt0, zt, params):
    
    # Unpack utility parameters and grids.
    alpha, beta, delta = params.alpha, params.beta, params.delta
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
# 3 - Perform value function iteration.
# -----------------------------------------------------------------------------------------------------
for iteration in range(1,number_of_iterations):
    
    # Loop over all possible starting states.
    for kt0 in enumerate(k_values):
        for zt in enumerate(z_values):
            
            # Solve Value Function and Policy Function and update values.
            V, g = Iterate_Value_Function(Value_Function[iteration-1,:,:], kt0, zt, params)
            
            Value_Function[iteration, kt0[0], zt[0]]  = V
            Policy_Function[iteration, kt0[0], zt[0]] = g
            
# -----------------------------------------------------------------------------------------------------
# 4 - Plot results.
# ----------------------------------------------------------------------------------------------------- 
# Plot Value Function for different starting states.
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for zt in enumerate(z_values):
    ax.plot(k_values, Value_Function[number_of_iterations-1,:,zt[0]])
ax.set(xlabel='k', ylabel='V(k,z)', title="Value Function")
ax.legend(np.flip(z_values).round(2), loc='right', title="z value")
plt.show()

# Plot final Policy Function for certain productivity values.
z_indices = [0,3,5,7,10]
fig, ax = plt.subplots()
for z_ind in z_indices:
    ax.plot(k_values, Policy_Function[number_of_iterations-1,:,z_ind])
ax.plot(k_values,k_values, '--', color='k', linewidth=0.8)
ax.set(xlabel='k', ylabel='g(k,z)', title="Policy Function")
ax.legend(np.flip(z_values[z_indices]).round(2), loc='right', title="z value")
plt.show()