

# assign initial values
α = 0.400
β = 0.987
δ = 0.012
iterations = 1000

# calculate steady-state level of capital
k_steady = ((1-β*(1-δ))/(α*β)) ^ (1/(α-1))                 # 100.44

# create deviations around steady state
k_values = range(0.98*k_steady, 1.02*k_steady, length=101)  # -2% to +2% of k_steady
k_values = round.(k_values, digits=2)                       # broadcast .

# create Bellman value function
# stores function value for each Kt value (column index), in each iteration (row index)
# example: Value_Function[50, 2] = Bellman value for second Kt from k_values on iteration 50
Value_Function = zeros(iterations, length(k_values))

# create optimal policy function
# for each input Kt (column index) store optimal Kt+1 index (row value), for each iteration (row index)
# example: Policy_Function[50, 2] = Index of optimal k_value for second k from k_values as input
#                                   for iteration 50
Policy_Function_Capital = zeros((iterations, length(k_values)))

# set Zt (productivity shock) equal to its expected value for now
zt = 1

#iteration=2
#kt0 = kt1 = (11, 100.44) 

# kt0[1] corresponds to the kt's index in k_value
# kt0[2] corresponds to the value of the kt

# execute value function iteration


for iteration in 2:iterations

    println(iteration) # track progress

    # loop though all possible starting states
    for kt0 in enumerate(k_values)

        v_max = nothing   # reset highest value function value

        # loop though all possible kt+1 choices
        for kt1 in enumerate(k_values)

            # calculate value of Value Function for given kt state and kt+1 choice
            new_value_function_value = log(zt*(kt0[2]^α) - kt1[2] + (1-δ)*kt0[2]) + β*Value_Function[iteration-1, kt1[1]]

            # check if there v_max value for this kt yet, set new one if not
            if v_max==nothing
                v_max = new_value_function_value
                # set max value for Value Function with this kt input
                Value_Function[iteration, kt0[1]] = v_max

            end

            # check if v_max from this kt+1 is a new maximum
            if new_value_function_value > v_max
                
                # update value of v_max for this kt state
                v_max = new_value_function_value
                # set max value for Value Function with this kt input
                Value_Function[iteration, kt0[1]] = v_max
                # store index of optimal kt+1 value in Policy Function
                Policy_Function_Capital[iteration, kt0[1]] = kt1[1]

            end
        end
    end
end
    

# plot results
using Plots, LaTeXStrings


# optimal policy function
optimal_k_index = Policy_Function_Capital[iterations,:] # get index of optimal k_values
optimal_k_values = k_values[round.(Int, optimal_k_index)] # get corresponding value

plot(k_values, optimal_k_values, label=L"g(k)", legend=:bottomright)
plot!(k_values,k_values, label=L"45^o Line", color=:black)

xlabel!("k")
ylabel!("g(k)")
title!("Policy Function")












