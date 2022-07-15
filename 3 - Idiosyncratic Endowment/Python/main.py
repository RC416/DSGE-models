"""
Idiosyncratic Entitlement model implemented in Python.

Steps:
    1 - 
    2 - 
    3 - 
    
"""

import numpy as np

from custom_functions import Solve_Value_Function, Get_Population_Distribution

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------

# Endowment parameters.
e_high = 1.0                                        # high entitlement
e_low  = 0.1                                        # low entitlement
e_grid = np.array([e_low, e_high])                  
number_of_e_values = 2
e_probs = np.array([[0.500, 0.500],[0.075, 0.925]]) # transition probabilities

# Utility parameters.
beta = 0.99322                                      # discounting factor (1/6th of year)
s = 1.5                                             # risk aversion coefficient

# Credit parameters.
a_high = 4                                          # upper credit limit
a_low  = -4                                          # lower credit limit
number_of_a_values = 100                            # credit grid size
a_grid = np.linspace(a_low, a_high, number_of_a_values)

# Store utility parameters and capital/productivity grids in a struct for passing to a function.
from typing import NamedTuple # could also use dataclass from dataclasses module

class Parameters(NamedTuple):
    s:       float
    beta:    float
    a_grid:  float
    e_grid:  float
    e_probs: float
    number_of_a_values: int
    number_of_e_values: int
    
params = Parameters(s, beta, a_grid, e_grid, e_probs, number_of_a_values, number_of_e_values)

# -----------------------------------------------------------------------------------------------------
# 2 - Solve model and get market bond price using binary search.
# -----------------------------------------------------------------------------------------------------

# Range of bond price values to search.
q_min = 0.8000
q_max = 1.2000
mcc   = np.Inf

# Iteration parameters.
dist = np.Inf
iteration_count = 0
max_iterations = 20
tolerance = 1e-4

# Solve for q.
while (dist > tolerance) & (iteration_count < max_iterations):
    
    # Get value of q from middle of range.
    q = (q_min + q_max)/2
    
    # Solve for Value Function and Policy Function
    Value_Function, Policy_Function, Policy_Function_Index = Solve_Value_Function(q, params)

    # Get Population Distribution.
    Population_Distribution = Get_Population_Distribution(Policy_Function_Index, params)

    # Check market clearing condition.
    mcc = (Policy_Function * Population_Distribution).sum()
    #if Policy_Function.sum() == 0: mcc = float('inf')

    # Update search parameters.
    dist = abs(mcc)
    iteration_count += 1
    
    # Update range of q according to the sign of mcc.
    if mcc > 0:  q_min = q
    if mcc < 0:  q_max = q
    
    # Print results.
    print(f"Iteration {iteration_count}: q={round(q,6)}, mcc={round(mcc,6)}")
    if iteration_count >= max_iterations: print(f"Warning: value function did not converge after {iteration_count} iterations")
    
# -----------------------------------------------------------------------------------------------------
# 3 - Plot results.
# ----------------------------------------------------------------------------------------------------- 

# Replicates Figure 1. of Huggett 1993.
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(a_grid, Policy_Function[:,1])
ax.plot(a_grid, Policy_Function[:,0])
ax.plot(a_grid, a_grid, '--', color='k', linewidth=0.8)
ax.set(xlabel="starting credit level", ylabel="optimal new credit level")
ax.legend(["high endowment", "low endowment", "45-degree line"], loc='lower right')
plt.show()

# Cumulative distribution of credit levels.
fig, ax = plt.subplots()
ax.plot(a_grid, Population_Distribution.sum(axis=1).cumsum())
ax.set(xlabel="starting credit level", title="Cumulative Distribution Function for Credit Level")
plt.show()












