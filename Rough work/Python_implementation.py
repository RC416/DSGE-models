# -*- coding: utf-8 -*-
"""
ECO2061
Assignment 2
Raymond Chiasson
1007278337

Question 3

Implementation of the MatLab template code
Uses the inputs from the assignment question with static Zt

"""

# load packages
import pandas as pd
import numpy as np
from math import log


# assign initial values
alpha = 0.400
beta  = 0.987
delta = 0.012
iterations = 1000


# create range of capital - must have positive log term of Bellman equation
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ** (1/(alpha-1))    # = 100.44 
k_values = [k_steady * x/10000 for x in range(9800, 10201, 4)]    #  -2% to +2% of k_steady 
k_values = [round(k,2) for k in k_values]

# create Bellman value function
# stores function value for each Kt value (column index), in each iteration (row index)
# example: Value_Function[50, 2] = Bellman value for second Kt from k_values on iteration 50
Value_Function = np.zeros((iterations, len(k_values)))


# create optimal policy function
# for each input Kt (column index) store optimal Kt+1 index (row value), for each iteration (row index)
# example: Policy_Function[50, 2] = Index of optimal k_value for second k from k_values as input
#                                   for iteration 50
Policy_Function_Capital = np.zeros((iterations, len(k_values)))


# set shock Zt to expected value for now
z=1


# execute value function iteration
for iteration in [x for x in range(1, iterations)]:

    # loop though all pairs of 'Kt' and 'Kt+1' 
    for kt0 in enumerate(k_values):
        
        # track the highest value of the Value Function each iteration
        v_max = None         # reset for every new value of Kt
        
        
        # loop though all possible Kt+1 choices for a given Kt
        for kt1 in enumerate(k_values):
    
            # calculate value of Value Function for given Kt state and Kt+1 choice
            New_Value_Function_Value = (log((z*(kt0[1] ** alpha)) - kt1[1] + (1-delta)*kt0[1])
                                        + beta*Value_Function[iteration-1, kt1[0]])
           
            
            
            # check if this Kt+1 value is the optimal choice for this Kt
            if (v_max is None) or (New_Value_Function_Value > v_max):
            
                # update value of v_max for this Kt state
                v_max = New_Value_Function_Value
            
                # set max value for Value Function with this Kt input
                Value_Function[iteration, kt0[0]] = v_max 
                
                # store index of optimal Kt+1 value in Policy Function
                Policy_Function_Capital[iteration, kt0[0]] = kt1[0]
        



# Policy_Function columns and values correspond to indeces of k_values
# convert these indeces to values of Kt (columns) and Kt+1 (values)

# use dataframe to better handle column names
Policy_Function_Capital_df = pd.DataFrame(data=Policy_Function_Capital,
                                          columns=[round(k, 2) for k in k_values])
# get k_value for each index
Policy_Function_Capital_df = Policy_Function_Capital_df.applymap(lambda x: round(k_values[int(x)],2))



'''
Notes on loop syntax:

enumerate() creates objects where
  kt0[1] is a value in the list k_values
  kt0[0] is the corresponding index.
      Example: kt1=(5, 100.44)

Varible types:

New_Value_Function_Value = float containing function output
Value_Function[iteration, Kt, Zt] = 3d array returning value for given state variables Kt, Zt
Policy_Function[iteration, Kt] = 2d array returning index of optimal Kt+1 for given Kt

'''


# Plot Value_Function iterations

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

for iteration in [x for x in range(1,iterations, int(iterations/10))]:
    ax.plot(k_values, Value_Function[iteration,])

ax.set(xlabel='kt', ylabel='V(kt)')

plt.show()



























