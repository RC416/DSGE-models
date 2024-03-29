"""
Idiosyncratic Endowment model implemented in Python.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Solve model and get market bond price using binary search.
        - Guess bond price
        - Solve model and get market bond price using binary search.
        - Get the distribution of credit and productivity levels
        - Check the market clearing condition and update bond price guess
    3 - Plot results.
"""

import numpy as np
from custom_functions import Solve_Value_Function, Get_Population_Distribution

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------

# Endowment parameters.
e_high = 1.0                                                # high endowment
e_low  = 0.1                                                # low endowment
e_grid = np.array([e_low, e_high])                          # endowment grid
number_of_e_values = 2
e_probs = np.array([[0.500, 0.500], [0.075, 0.925]])         # transition probabilities

# Utility parameters.
sigma = 1.5                                                 # risk aversion coefficient
beta = 0.99322                                              # discount factor

# Credit parameters.
a_high = 4                                                  # upper credit limit
a_low  = -2                                                 # lower credit limit / borrowing constraint
number_of_a_values = 100                                    # credit grid size
a_grid = np.linspace(a_low, a_high, number_of_a_values)     # credit grid

# Store parameters in a struct for passing to a function.
from typing import NamedTuple # could also use dataclass from dataclasses module

class Parameters(NamedTuple):
    sigma:   float
    beta:    float
    a_grid:  float
    e_grid:  float
    e_probs: float
    number_of_a_values: int
    number_of_e_values: int
    
params = Parameters(sigma, beta, a_grid, e_grid, e_probs, number_of_a_values, number_of_e_values)

# -----------------------------------------------------------------------------------------------------
# 2 - Solve model and get market bond price using binary search.
# -----------------------------------------------------------------------------------------------------

# Range of bond price values to search.
q_min = 0.985
q_max = 1.100

# Optional: floor for q_min such that those at credit limit can still afford positive consumption.
q_min = (a_low + e_low) / a_low

# Placeholder for the market clearing condition.
mcc = np.Inf

# Iteration parameters.
dist = np.Inf
iteration_count = 0
max_iterations = 20
tolerance = 1e-3

# Solve for market price q.
while (dist > tolerance) & (iteration_count < max_iterations):
    
    # Get value of q from middle of range.
    q = (q_min + q_max)/2
    
    # Solve for the Value Function and Policy Function.
    Value_Function, Policy_Function, Policy_Function_Index = Solve_Value_Function(q, params)

    # Get the Population Distribution.
    Population_Distribution = Get_Population_Distribution(Policy_Function_Index, params)

    # Check the market clearing condition.
    mcc = (Policy_Function * Population_Distribution).sum()

    # Update search parameters.
    dist = abs(mcc)
    iteration_count += 1
    
    # Update the range of q according to the sign of the market clearing condition.
    if mcc > 0:  q_min = q
    if mcc < 0:  q_max = q
    
    # Print results.
    print(f"Iteration {iteration_count}: q={round(q,6)}, mcc={round(mcc,6)}")
    if iteration_count >= max_iterations: print(f"Warning: search for q did not converge after {iteration_count} iterations")
    
# -----------------------------------------------------------------------------------------------------
# 3 - Plot results.
# ----------------------------------------------------------------------------------------------------- 

# Policy functions (Figure 1. from Huggett 1993).
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(a_grid, Policy_Function[:, 1])
ax.plot(a_grid, Policy_Function[:, 0])
ax.plot(a_grid, a_grid, '--', color='k', linewidth=0.8)
ax.set(xlabel="starting credit level", ylabel="optimal new credit level", title="Policy Function")
ax.legend(["high endowment", "low endowment", "45-degree line"], loc='lower right')
plt.show()

# Distribution of credit levels (Figure 2. from Huggett 1993).
fig, ax = plt.subplots()
ax.plot(a_grid, Population_Distribution[:,1].cumsum())
ax.plot(a_grid, Population_Distribution[:,0].cumsum())
ax.set(xlim=[-2, 1], ylim=[0, 1], xlabel="starting credit level", title="Cumulative Distribution Function for Credit Level")
ax.legend(["high endowment", "low endowment"], loc='upper left')
plt.show()