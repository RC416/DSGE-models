#=
Stochastic Growth model implemented in Python.
=#

# Assign parameter values.
α = 0.400;
β = 0.987;
δ = 0.012;
number_of_iterations = 1000;

# Calculate the steady-state level of capital.
k_steady = ((1-β*(1-δ))/(α*β)) ^ (1/(α-1))

# Create a range of capital values around steady-state (+/- 50%).
number_of_k_values = 201;
k_low_pct = 0.98;
k_high_pct = 1.02;
k_values = range(k_low_pct*k_steady, k_high_pct*k_steady, length=number_of_k_values);

# Get productivity levels and transition probabilities.
using DelimitedFiles
cd("C:\\Users\\Ray\\Documents\\GitHub\\DSGE-models\\2 - Stochastic Growth\\Julia implementation")
z_probs = readdlm("Inputs\\z_probs.csv", ',', Float64)
z_values = readdlm("Inputs\\z_values.csv", ',', Float64)[:,1]
number_of_z_values = size(z_values,1);

# Initialize Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);

# Perform value function iteration.
for iteration in 2:number_of_iterations

    # Loop over all possible starting states.
    for kt0 in eachindex(k_values)
        for zt in eachindex(z_values)

            # Solve Value Function and Policy Function and update values.
            V, g = Iterate_Value_Function_v2(Value_Function[iteration-1,:,:], kt0, zt);

            Value_Function[iteration, kt0, zt] = V;
            Policy_Function[iteration, kt0, zt] = g;
        end
    end
end


# Plot Value Function for different starting states.
using Plots
Figure1 = plot(legendtitle="z value")
for zt_index in eachindex(z_values)
    plot!(k_values, Value_Function[number_of_iterations,:,zt_index], label=round(z_values[zt_index],digits=2))
end
xlabel!("k")
ylabel!("V(k,z)")
title!("Value Function")
display(Figure1)

# Plot final Policy Function for certain productivity values.
z_indices = [1,4,6,8,11]
Figure2 = plot(k_values, Policy_Function[number_of_iterations,:,z_indices], label=round.(z_values[z_indices], digits=2)', legend=:topleft, legendtitle="z value")
plot!(k_values,k_values, label="45-degree line", color=:black)
xlabel!("k")
ylabel!("g(k,z)")
title!("Policy Function")
display(Figure2)
















#=
Function to perform an iteration of value function iteration.
Two different versions with increasing levels of performance.
Each version has idential inputs and outputs.

Input:
    - Current value function
    - Starting capital and productivity level 
    
Output:
    - Next value function value
    - Optimal choice of next period capital
    
Version
    v1: using only base functions + for-loops
    v2: using broadcast/vectorized calculation instead of for-loop
=#

using LinearAlgebra

#=
Version 1 - using only base functions + for-loops.
=#

function Iterate_Value_Function_v1(Previous_Value_Function, kt0_index, zt_index)

    # Get capital and productivity values from index.
    kt0 = k_values[kt0_index];
    zt = z_values[zt_index];

    v_max = -Inf;
    kt1_optimal = 0.0;

    # Search over possible next period capital choices.
    for kt1_index in eachindex(k_values)

        # Get capital value from index.
        kt1 = k_values[kt1_index];

        # Calculate value function for choice of next period capital.
        new_v_max = log(zt*(kt0^α) - kt1 + (1-δ)*kt0) + β*dot(Previous_Value_Function[kt1_index,:],z_probs[zt_index,:]);

        # Check if this capital choice gives highest Value Function value.
        if new_v_max > v_max
        
            # Update candidate values.
            v_max = new_v_max;
            kt1_optimal = kt1;
        end
    end

    return v_max, kt1_optimal;
end


#=
Version 2 - using broadcast/vectorized calculation instead of for-loop.
=#

function Iterate_Value_Function_v2(Previous_Value_Function, kt0_index, zt_index)

    # Get capital and productivity values from index.
    kt0 = k_values[kt0_index];
    zt = z_values[zt_index];

    # Calculate array of value function values for all k_values.
    V_max_values = log.(zt*(kt0^α) + (1-δ)*kt0 .- k_values) + β*(Previous_Value_Function*z_probs[zt_index,:]);

    # Get index for the optimal capital choice.
    kt1_index_optimal = argmax(V_max_values);

    # Get Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_index_optimal];
    v_max = V_max_values[kt1_index_optimal];

    return v_max, kt1_optimal;
end


# Benchmarking
using BenchmarkTools
@btime Iterate_Value_Function_v1(Value_Function[number_of_iterations,:,:], 1, 1) # 180 μs
@btime Iterate_Value_Function_v2(Value_Function[number_of_iterations,:,:], 1, 1) # 8.2 μs