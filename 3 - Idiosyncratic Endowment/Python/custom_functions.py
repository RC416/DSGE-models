# -*- coding: utf-8 -*-
"""
Three functions with multiple versions:
    1. Solve Household Problem
        - v1, v2
    2. Solve Value Function
        - Uses Solve Household Problem
    3. Get Population Distribution

-----------------------------------------------------------------------------------------------------
Function 1: Solve Household Problem
    
Function to perform an iteration of value function iteration.
Four different versions with increasing levels of performance.
Each version has identical inputs and outputs.

Input:
    - Current Value Function
    - Market price for bond
    - Starting wealth and endowment level 
    
Output:
    - Next value function value
    - Optimal choice of savings/borrowing for next period
    
Version
    v1: using only base functions + for-loops
    v2: using Numpy arrays / broadcast instead of for-loop
    Using numba @jit precopiling is optional and compatible with both.
    Must also use numba @jit with Solve Value Function.

-----------------------------------------------------------------------------------------------------
Function 2: Solve Value Function

Solves for the Value Function and Policy Function using value function iteration.
Applies Solve Household Problem to all possible starting states for each iteration.
Search parameters (max iterations, tolerance, etc.) are defined in the function.

Input:
    - Market price for bond
    
Output:
    - Value Function
    - Policy Function
    - Policy Function (with index values)
        
-----------------------------------------------------------------------------------------------------
Function 3: Get Population Distribution. 

Solves for the steady state distribution over credit-endowment states given
a policy function (with index values).
Search parameters (max iterations, tolerance, etc.) are defined in the function.

Input:
    - Policy Function (with index values)

Output:
    - Steady-state population distribution

"""

import numpy as np

# -----------------------------------------------------------------------------------------------------
# Function 1 - Version 1 - Solve Household Problem using only base functions + for-loops.
# -----------------------------------------------------------------------------------------------------
from numba import jit

@jit(nopython=True)
def Solve_HH_Problem_v1(q, a_start_index, e_start_index, Value_Function, params):

    # Unpack utility parameters and grids.
    s, beta = params.s, params.beta
    a_grid  = params.a_grid
    e_grid  = params.e_grid
    e_probs = params.e_probs
    
    # Get value of state variables.
    a_start = a_grid[a_start_index]
    e_start = e_grid[e_start_index]
    
    # Variables to store candidate optimal values for Value Function and Policy Function.
    v_max = -np.Inf
    a_next_optimal = a_grid[0]
    a_next_index = 0
    
    # Search over possible next period borrowing choices.
    for a_next in enumerate(a_grid):
            
        # Check budget constraint: if consumption is negative, skip this value.
        if (a_start + e_start - q*a_next[1] <= 0):
            break
        
        # Calculate value function for given choice of next period capital.
        else:
            new_v_max = ((a_start + e_start - q*a_next[1]) ** (1 - s))/(1-s) + beta*np.dot(Value_Function[a_next[0],:], e_probs[e_start_index,:])
            
        # Check if this capital choice gives highest Value Function value.
        if new_v_max > v_max:
        
            # Update candidate values.
            v_max = new_v_max
            a_next_optimal = a_next[1]
            a_next_index   = a_next[0]
    
    return v_max, a_next_optimal, a_next_index

# -----------------------------------------------------------------------------------------------------
# Function 1 - Version 2 - Solve Household Problem using vectorized calculation
# -----------------------------------------------------------------------------------------------------
from numba import jit

@jit(nopython=True)
def Solve_HH_Problem_v2(q, a_start_index, e_start_index, Value_Function, params):

    # Unpack utility parameters and grids.
    s, beta = params.s, params.beta
    a_grid  = params.a_grid
    e_grid  = params.e_grid
    e_probs = params.e_probs
    
    # Get value of state variables.
    a_start = a_grid[a_start_index]
    e_start = e_grid[e_start_index]
       
    # Vector of consumption values dictated by credit selection.
    Consumption = a_start + e_start - q*a_grid
    valid_indices = (Consumption > 0)
    
    # Calculate value function values.
    V_max_values = np.power(Consumption[valid_indices],(1-s)) / (1-s) + beta*np.dot(Value_Function[(valid_indices),:], e_probs[e_start_index,:])

    # Get index of optimal credit choice within the vector of valid indices.
    optimal_subindex = np.argmax(V_max_values)
    
    # Get values and original index of optimal values.
    a_next_optimal = (a_grid[valid_indices])[optimal_subindex]
    v_max = V_max_values[optimal_subindex]
    a_next_index = np.where(a_grid == a_next_optimal)[0][0]
    
    return v_max, a_next_optimal, a_next_index

# -----------------------------------------------------------------------------------------------------
# Function 2 - Solve for Value Function and Policy Function.
# -----------------------------------------------------------------------------------------------------
from numba import jit

@jit(nopython=True)
def Solve_Value_Function(q, params):
       
    # Unpack relevant parameters.
    number_of_a_values = params.number_of_a_values
    number_of_e_values = params.number_of_e_values
    
    # Arrays to hold 2 value function iterations.
    Value_Function = np.zeros((number_of_a_values, number_of_e_values))
    Policy_Function = np.zeros((number_of_a_values, number_of_e_values))
    Policy_Function_Index = np.zeros((number_of_a_values, number_of_e_values), dtype=np.int64)
    
    Value_Function_New = np.zeros((number_of_a_values, number_of_e_values))
    Policy_Function_New = np.zeros((number_of_a_values, number_of_e_values))
    Policy_Function_Index_New = np.zeros((number_of_a_values, number_of_e_values), dtype=np.int64)
    
    # Iteration parameters.
    dist = np.Inf
    iteration_count = 0
    max_iterations = 5000
    tolerance = 1e-6
    
    # Solve for Value Function and Policy Function.
    while (dist > tolerance) & (iteration_count < max_iterations):
              
        # Loop over all possible starting states.
        for a_start_index in range(number_of_a_values):
            for e_start_index in range(number_of_e_values):
                            
                # Solve Value Function and Policy Function and update values.
                V, g, g_index = Solve_HH_Problem_v1(q, a_start_index, e_start_index, Value_Function, params) 
                #V, g, g_index = Solve_HH_Problem_v2(q, a_start_index, e_start_index, Value_Function, params)


                Value_Function_New[a_start_index, e_start_index] = V
                Policy_Function_New[a_start_index, e_start_index] = g
                Policy_Function_Index_New[a_start_index, e_start_index] = g_index
        
        # Update search parameters.
        dist = np.abs(Value_Function - Value_Function_New).max()
        iteration_count += 1
        
        # Update Value Function.
        Value_Function = Value_Function_New.copy()
        Policy_Function = Policy_Function_New.copy()
        Policy_Function_Index = Policy_Function_Index_New.copy()
        
        # Print warning if convergence is not achieved.
        if iteration_count >= max_iterations:
            print(f"Warning: value function did not converge after {iteration_count} iterations")
    
    return Value_Function, Policy_Function, Policy_Function_Index

# -----------------------------------------------------------------------------------------------------
# Function 3 - Version 1 - Get population distribution from policy function
# -----------------------------------------------------------------------------------------------------
from numba import jit

@jit(nopython=True)
def Get_Population_Distribution(Policy_Function_Index, params):
      
    # Unpack key varaibles.
    number_of_a_values = params.number_of_a_values
    number_of_e_values = params.number_of_e_values
    e_probs = params.e_probs
    
    # Optional: create local copy of e_probs that is contiguous on column slices (numerically identical).
    #e_probs = e_probs.copy(order='F') 

    # Arrays to store 2 iterations of finding population distribution.    
    Population_Distribution = np.ones((number_of_a_values, number_of_e_values)) / (number_of_a_values * number_of_e_values)
    New_Distribution = Population_Distribution.copy()
    

    # Iteration parameters.
    dist = np.Inf
    iteration_count = 0
    max_iterations = 5000
    tolerance = 1e-10
    
    # Solve for population distribution.
    while (dist > tolerance) & (iteration_count < max_iterations):
        
        # Get "inflow" to each credit-endowment state in next period
        for a_index in range(number_of_a_values):
            for e_index in range(number_of_e_values):
                
                # Get indices of states that flow into given state.               
                inflow_indices = np.zeros((number_of_a_values, number_of_e_values))
                for a_index_inflow in range(number_of_a_values):
                    for e_index_inflow in range(number_of_e_values):
                        if Policy_Function_Index[a_index_inflow, e_index_inflow] == a_index:
                            inflow_indices[a_index_inflow, e_index_inflow] = 1

                # Alternative if not using numba precompiling. Approximately the same speed as numba.
                #inflow_indices = (Policy_Function_Index == a_index).astype(float)                
                
                # Sum distribution-weighted inflow into given state.
                inflow = (np.multiply(Population_Distribution, e_probs[:,e_index]) * inflow_indices).sum()    
                
                New_Distribution[a_index, e_index] = inflow

        # Update search parameters.
        dist = np.abs(Population_Distribution - New_Distribution).max()
        iteration_count += 1
        
        # Update population distribution.
        Population_Distribution = New_Distribution.copy()
                
        # Print warning if convergence is not achieved.
        if iteration_count >= max_iterations:
            print(f"Warning: population distribution did not converge after {iteration_count} iterations")
            
    return Population_Distribution