#=
Stochastic Growth model implemented in Julia.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Create function to solve Value Function for given starting states.
    3 - Perform value function iteration.
    4 - Plot results.
=#

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------
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
k_values = collect(range(k_low_pct*k_steady, k_high_pct*k_steady, length=number_of_k_values));

# Get productivity levels and transition probabilities.
using DelimitedFiles
cd("C:\\Users\\Ray\\Documents\\GitHub\\DSGE-models\\2 - Stochastic Growth\\Julia implementation")
z_probs = readdlm("Inputs\\z_probs.csv", ',', Float64);
z_values = readdlm("Inputs\\z_values.csv", ',', Float64)[:,1];
number_of_z_values = size(z_values,1);

# Initialize Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);

# Store utility parameters and capital/productivity grids in a struct for passing to a function.
struct Parameters
    α::Float64
    β::Float64
    δ::Float64
    k_values::Array{Float64, 1} # == Vector{Float64}
    z_values::Array{Float64, 1} # == Vector{Float64}
    z_probs::Array{Float64, 2}  # == Matrix{Float64} 
end

params = Parameters(α,β,δ,k_values,z_values,z_probs);

# -----------------------------------------------------------------------------------------------------
# 2 - Create function to solve Value Function for given starting states.
# -----------------------------------------------------------------------------------------------------
function Iterate_Value_Function(Previous_Value_Function, kt0_index, zt_index, params)

    # Unpack utility parameters and grids.
    α,β,δ = params.α, params.β, params.δ;
    k_values = params.k_values;
    z_values = params.z_values;
    z_probs = params.z_probs;

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

# Alternative: import other versions from module.
#include("iteration_functions.jl")
#using .iteration_functions

# -----------------------------------------------------------------------------------------------------
# 3 - Perform value function iteration.
# -----------------------------------------------------------------------------------------------------
for iteration in 2:number_of_iterations

    # Loop over all possible starting states.
    for kt0_index in eachindex(k_values), zt_index in eachindex(z_values)

        # Solve Value Function and Policy Function and update values.
        V, g = Iterate_Value_Function(Value_Function[iteration-1,:,:], kt0_index, zt_index, params);

        Value_Function[iteration, kt0_index, zt_index] = V;
        Policy_Function[iteration, kt0_index, zt_index] = g;
    end
end

# -----------------------------------------------------------------------------------------------------
# 4 - Plot results.
# -----------------------------------------------------------------------------------------------------
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
