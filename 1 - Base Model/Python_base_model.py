# -*- coding: utf-8 -*-
"""
Base Model implementation in Python
"""

# Load packages.
import numpy as np
from math import log

# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 1.000
number_of_iterations = 1000

# Calculate the steady-state level of capital.
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ** (1/(alpha-1))    

# Create a range of capital values around steady-state (+/- 50%).
number_of_k_values = 201
k_low_pct = 0.50
k_high_pct = 1.50
k_values = np.linspace(k_low_pct*k_steady, k_high_pct*k_steady, num=number_of_k_values)
#k_values = np.round(k_values, decimals=2, out=None)

# Initialize Value Function and Policy Function (as arrays).
Value_Function = np.zeros((number_of_iterations, number_of_k_values))
Policy_Function = np.zeros((number_of_iterations, number_of_k_values))


# Perform value function iteration.
for iteration in range(1,number_of_iterations):
    
    for kt0 in enumerate(k_values):         # for each level of starting capital...
        
        v_max = None # placeholder for candidate value function value
        
        for kt1 in enumerate(k_values):     # ...check all next period capital choices 
    
            # Calculate value of Value Function for given starting capital and next period capital choice.
            New_Value_Function_Value = (log((kt0[1] ** alpha) - kt1[1] + (1-delta)*kt0[1])
                                        + beta*Value_Function[iteration-1, kt1[0]])
           
            # Check if this capital choice gives highest Value Function value.
            if (v_max is None) or (New_Value_Function_Value > v_max):
            
                # Update value of v_max, Value Function, and Policy Function
                v_max = New_Value_Function_Value
                Value_Function[iteration, kt0[0]] = v_max       # store value
                Policy_Function[iteration, kt0[0]] = kt1[0]     # store index


# Plot Value Function for various iterations.
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for iteration in range(0,number_of_iterations, int(number_of_iterations/10)):
    ax.plot(k_values, Value_Function[iteration,])

ax.set(xlabel='kt', ylabel='V(kt)', title="Value Function")
ax.legend(range(0,number_of_iterations, int(number_of_iterations/10)), loc='right')
plt.show()


# Plot final Policy Function.
fig, ax = plt.subplots()
ax.plot(k_values, k_values[Policy_Function[-1,].astype(int)])
ax.plot([min(k_values),max(k_values)],[min(k_values),max(k_values)], '--', color='k', linewidth=0.8)
ax.set(xlabel='k', ylabel='g(k)', title="Policy Function")
ax.legend(["g(k)", "45-degree line"])
plt.show()


'''
Notes on syntax:
    
The Value Function and Policy Function are stored as arrays where 
the indices correspond to the iteration number and the starting 
level of capital (kt), resepectively. 

enumerate() creates tuples for each element in a list where
  kt0[1] is a value in the list k_values
  kt0[0] is the corresponding index.
      Example: kt1=(5, 100.44)

'''