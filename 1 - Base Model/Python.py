# Base Model implementation in Python.

# Load packages.
import numpy as np
from math import log

# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 1.000
number_of_iterations = 1000

# Calculate the steady-state level of capital.
k_steady = ((1 - beta * (1 - delta)) / (alpha * beta)) ** (1 / (alpha - 1))    

# Create a grid of capital values around the steady-state (+/- 50%).
number_of_k_values = 201
k_low_pct = 0.50
k_high_pct = 1.50
k_values = np.linspace(k_low_pct * k_steady, k_high_pct * k_steady, num=number_of_k_values)

# Initialize the Value Function and Policy Function (as arrays).
Value_Function = np.zeros((number_of_iterations, number_of_k_values))
Policy_Function = np.zeros((number_of_iterations, number_of_k_values))

# Solve the household's problem for each possible starting state.
for iteration in range(1, number_of_iterations):
    
    for kt0 in enumerate(k_values):
        
        # Variables to store candidate optimal values.
        v_max = float('-inf')
        kt1_optimal = 0.0
        
        # Check all possible next period capital choices.
        for kt1 in enumerate(k_values): 
    
            # Calculate the Value Function for given starting capital and next period capital choice.
            new_value_function_value = (log((kt0[1] ** alpha) - kt1[1] + (1 - delta) * kt0[1])
                                        + beta * Value_Function[iteration - 1, kt1[0]])
           
            # Check if this capital choice gives the highest Value Function value.
            if new_value_function_value > v_max:
            
                # Update candidate values.
                v_max = new_value_function_value
                kt1_optimal = kt1[1]
                
        # Update the Value Function and Policy Function with optimal values. 
        Value_Function[iteration, kt0[0]] = v_max
        Policy_Function[iteration, kt0[0]] = kt1_optimal

# Plot various iterations of the Value Function.
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for iteration in range(0,number_of_iterations, int(number_of_iterations / 10)):
    ax.plot(k_values, Value_Function[iteration, ])
ax.set(xlabel='k', ylabel='V(k)', title="Value Function")
ax.legend(range(1, number_of_iterations + 1, int(number_of_iterations / 10)), loc='right')
plt.show()

# Plot the final Policy Function.
fig, ax = plt.subplots()
ax.plot(k_values, Policy_Function[-1,])
ax.plot(k_values, k_values, '--', color='k', linewidth=0.8)
ax.set(xlabel='k', ylabel='g(k)', title="Policy Function")
ax.legend(["g(k)", "45-degree line"])
plt.show()
