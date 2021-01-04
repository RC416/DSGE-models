# -*- coding: utf-8 -*-
"""
ECO2061
Assignment 2
Raymond Chiasson
1007278337

Solution to Question 3

Solve the value function using zt shocks.

"""

# load packages
import pandas as pd
import numpy as np
from math import log


'''
1 - Load data and set up variables
'''

# assign initial values
alpha = 0.400
beta  = 0.987
delta = 0.012
rho   = 0.950
sigma = 0.007

iterations = 100
final_iteration = iterations - 1


# create range of capital
# must have positive log term of Bellman equation
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ** (1/(alpha-1))   # = 100.44 
k_values = [k_steady * x/10000 for x in range(9800, 10201, 2)]    #  -2% to +2% of k_steady 


# get values of Z and transition probabilities from Excel sheets
z_values = pd.read_excel('Z.xlsx', header=None, index_col=0)
z_probs  = pd.read_excel('Zprob.xlsx', header=0, index_col=0)     # from "row state" to "col state"


# create Bellman value function
# stores a value for each iteration, kt, zt point
# example: Value_Function[50, 2, 5] = value for iteration 50, second kt, fifth zt
Value_Function = np.zeros((iterations, len(k_values), len(z_values)))


# create optimal policy function
# for each input Kt (column index) store optimal Kt+1 index (row value), for each iteration (row index)
# example: Policy_Function[50, 2] = Index of optimal k_value for second k from k_values as input
#                                   for iteration 50
Policy_Function_Capital = np.zeros((iterations, len(k_values), len(z_values)))



'''
2 - Create function to calculate/retrieve value function value
'''

def get_value_function_value(zt0, kt0, kt1):
    
    # value fuction given by:
    # V(kt, zt) = log(ct) + Beta * E[V(kt+1, zt+1)]
    # where log(ct) = current period,
    #  and E[V(kt+1, zt+1)] = future periods

    # get the value of the current period term
    current_period = log((zt0[1]*(kt0[1] ** alpha)) - kt1[1] + (1-delta)*kt0[1])
    
    
    # calculate the expected value for future periods given zt1 transitions
    future_periods = 0  
    
    # loop though all possible shocks next period
    for zt1 in enumerate(z_values[1]):
        
        # get probability of given zt1 shock, given zt0 stating point
        probability = z_probs.iloc[zt0[0], zt1[0]]
        # get value function value for given zt1 shock
        VF_value = Value_Function[iteration-1, kt1[0], zt1[0]]
        
        # add expected value to future periods value
        future_periods += probability * VF_value
        
        
    return current_period + beta*future_periods

# zt0, kt0, and kt1 come from the loop below. See notes below loop.
# kt[0] = index
# kt[1] = value
# same for kt0, zt0



''' 
3 - Solve the value function
'''

# let iteration 0 be the all 0s guess
for iteration in [x for x in range(1, iterations)]:

    print(iteration) # track progress

    # loop though all starting shocks zt0
    for zt0 in enumerate(z_values[1]):

        # loop though all pairs of Kt and Kt+1
        for kt0 in enumerate(k_values):
            
            # track the highest value of the Value Function each iteration
            v_max = None
            
            
            # loop though all possible Kt+1 choices for a given Kt
            for kt1 in enumerate(k_values):
        
                # calculate value of Value Function for given Kt/Zt state and Kt+1 choice          
                New_Value_Function_Value = get_value_function_value(zt0, kt0, kt1)
                
                # check if this Kt+1 value is the optimal choice for this Kt
                if (v_max is None) or (New_Value_Function_Value > v_max):
                
                    # update value of v_max for this (Kt,Zt) state
                    v_max = New_Value_Function_Value
                
                    # set max value for Value Function with this (Kt, Zt) input
                    Value_Function[iteration, kt0[0], zt0[0]] = v_max 
                    
                    # store index of optimal Kt+1 value in Policy Function
                    Policy_Function_Capital[iteration, kt0[0], zt0[0]] = kt1[0]
            
'''
Notes on loop syntax:

enumerate() creates objects where
  kt0[1] is a value in the list k_values
  kt0[0] is the corresponding index.
      Example: kt1=(5, 100.44)

Varible types:

New_Value_Function_Value = float containing function output
Value_Function[iteration, Kt, Zt] = 3d array returning value for given state variables Kt, Zt
Policy_Function[iteration, Kt, Zt] = 3d array returning index of optimal Kt+1 for given Kt, Zt

'''



'''
4 - Optimal consumption policy function
'''

# get optimal Kt+1 from last iteration of policy function for capital
# here row=Kt, column=Zt, cell/value=Kt+1 (all indeces not actual value)
# create copy to avoid editing the original array
Policy_Function_Consumption = np.copy(Policy_Function_Capital)[final_iteration,:,:]


# function to calculate the resulting consumption choice from optimal Kt+1 allocation
def calculate_consumption(Kt0, Zt0, Kt1):
    # from the budget constraint: ct = f(zt, kt) - kt+1 + (1-delta)kt
    consumption = Zt0*(Kt0 ** alpha) - Kt1 + (1-delta)*Kt0    
    return consumption


# calculate optimal consumption choice for all sets of Kt, Zt
for kt0 in enumerate(k_values):
    
    for zt0 in enumerate(z_values[1]):
        
        # get index of kt+1 from policy function
        kt1_i = int(Policy_Function_Capital[final_iteration, kt0[0], zt0[0]])
        
        # calculate optimal consumption given Kt, Zt, Kt+1
        optimal_consumption = calculate_consumption(kt0[1], zt0[1], k_values[kt1_i])
        
        Policy_Function_Consumption[kt0[0], zt0[0]] = optimal_consumption


# convert to dataframe with row/column labels
Policy_Function_Consumption_df = pd.DataFrame(data=Policy_Function_Consumption,
                                              columns=[round(z, 2) for z in z_values[1]],
                                              index=[round(k, 2) for k in k_values])

# save to csv
Policy_Function_Consumption_df.to_csv('consumption_policy_function.csv')


'''
5 - Optimal capital policy function
'''
# convert to dataframe with row/column labels. Use final iteration.
Policy_Function_Capital_df = pd.DataFrame(data=Policy_Function_Capital[final_iteration,:,:],
                                          columns=[round(z, 2) for z in z_values[1]],
                                          index=[round(k, 2) for k in k_values])
# get k_value for each index
Policy_Function_Capital_df = Policy_Function_Capital_df.applymap(lambda x: round(k_values[int(x)],2))

# save to csv
Policy_Function_Capital_df.to_csv('capital_policy_function.csv')



'''
6 - Plot Value Functions
'''

import matplotlib.pyplot as plt

# create ploting objects
fig, ax = plt.subplots()

# plot value function curve for each Zt
for Zt0 in enumerate(z_values[1]):
    ax.plot(k_values, Value_Function[final_iteration,:,Zt0[0]])

# label axes
ax.set(xlabel='$k_t$', ylabel='V($k_t$, $z_t$)',
       title='Value Function for all $k_t$ and $z_t$')

# create legend
lgd = ax.legend(z_values.index,
                loc='right', bbox_to_anchor=(1.25,0.6), # position legend
                labelspacing=-2.6, frameon=False) # to reverse legened

# display and save figure
plt.show()
fig.savefig('Question 3 - value function plot v2.png', format='png', dpi=1200,
            bbox_extra_artists=(lgd,), bbox_inches='tight') # to get legend in the frame



'''
7 - Plot Policy Functions
'''

import matplotlib.pyplot as plt

# create ploting objects
fig, ax = plt.subplots()

# plot value function for each Zt
for Zt0 in enumerate(z_values[1]):
    ax.plot(k_values, Policy_Function_Capital_df.loc[:,round(Zt0[1],2)])

# label axes
ax.set(xlabel='$k_t$', ylabel='g($k_t$, $z_t$)',
       title='Capital Policy Function for all $k_t$ and $z_t$')

# create legend
lgd = ax.legend(z_values.index,
                loc='right', bbox_to_anchor=(1.0,0.35),   # position legend
                labelspacing=-2.6, frameon=False,        # to reverse legened
                fontsize=7) 

# add a 45 degree line
ax.plot([min(k_values),max(k_values)], [min(k_values),max(k_values)], 'g:', color='k', linewidth=1.3)

# display and save figure
plt.show()
fig.savefig('Question 3 - policy function plot v2.png', format='png', dpi=1200,
            bbox_extra_artists=(lgd,), bbox_inches='tight') # to get legend in the frame


























