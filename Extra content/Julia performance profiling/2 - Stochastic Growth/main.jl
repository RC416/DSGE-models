#=
Stochastic Growth model implemented in Julia.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Create a function to solve the household's problem for a given starting state.
    3 - Perform value function iteration.
    4 - Plot results.
=#

# Import modules and declare struct.
using DelimitedFiles
using Plots
cd("C:\\Users\\Ray\\Documents\\GitHub\\DSGE-models\\Extra content\\Julia performance profiling\\2 - Stochastic Growth")
include("custom_functions.jl")
using .custom_functions

# Store utility parameters and capital/productivity grids in a struct for passing to a function.
struct Parameters
    α::Float64
    β::Float64
    δ::Float64
    k_values::Array{Float64, 1} # == Vector{Float64}
    z_values::Array{Float64, 1} # == Vector{Float64}
    z_probs::Array{Float64, 2}  # == Matrix{Float64} 
end


function main()
# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------
# Assign parameter values.
α = 0.400;
β = 0.987;
δ = 0.012;
number_of_iterations = 1000;

# Calculate the steady-state level of capital.
k_steady = ((1 - β * (1 - δ)) / (α * β)) ^ (1 / (α - 1))

# Create a range of capital values around the steady-state (+/- 50%).
number_of_k_values = 201;
k_low_pct = 0.98;
k_high_pct = 1.02;
k_values = collect(range(k_low_pct * k_steady, k_high_pct * k_steady, length=number_of_k_values));

# Get productivity levels and transition probabilities.
z_probs = readdlm("Inputs\\z_probs.csv", ',', Float64);
z_values = readdlm("Inputs\\z_values.csv", ',', Float64)[:, 1];
number_of_z_values = size(z_values, 1);

# Initialize the Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values, number_of_z_values);

params = Parameters(α, β, δ, k_values, z_values, z_probs);

# -----------------------------------------------------------------------------------------------------
# 3 - Perform value function iteration.
# -----------------------------------------------------------------------------------------------------
for iteration in 2:number_of_iterations

    # Loop over all possible starting states.
    #for kt0_index in eachindex(k_values), zt_index in eachindex(z_values) # single-thread implementation
    Threads.@threads for zt_index in eachindex(z_values)                   # multi-thread implementation
        for kt0_index in eachindex(k_values)                               # multi-thread implementation

        # Solve the Value Function and Policy Function and update values.
        @views V, g = Solve_HH_Problem_v1(Value_Function[iteration-1, :, :], kt0_index, zt_index, params);
        #@views V, g = Solve_HH_Problem_v2(Value_Function[iteration-1, :, :], kt0_index, zt_index, params);

        Value_Function[iteration, kt0_index, zt_index] = V;
        Policy_Function[iteration, kt0_index, zt_index] = g;
        end
    end
end

# -----------------------------------------------------------------------------------------------------
# 4 - Plot results.
# -----------------------------------------------------------------------------------------------------
# Plot the Value Function for different starting states.
Figure1 = plot(legendtitle="z value")
for zt_index in eachindex(z_values)
    plot!(k_values, Value_Function[number_of_iterations, :, zt_index], label=round(z_values[zt_index], digits=2))
end
xlabel!("k")
ylabel!("V(k,z)")
title!("Value Function")
display(Figure1)

# Plot the final Policy Function for certain productivity values.
z_indices = [1,4,6,8,11]
Figure2 = plot(k_values, Policy_Function[number_of_iterations, :, z_indices], label=round.(z_values[z_indices], digits=2)', legend=:topleft, legendtitle="z value")
plot!(k_values,k_values, linestyle=:dash, label="45-degree line", color=:black)
xlabel!("k")
ylabel!("g(k,z)")
title!("Policy Function")
display(Figure2)

end

# Benchmarking and profiling.
using Profile
using BenchmarkTools

@btime main()
@profview main()
# https://www.julia-vscode.org/docs/dev/userguide/profiler/