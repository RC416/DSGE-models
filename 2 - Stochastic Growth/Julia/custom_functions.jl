#=
Function to solve the household's problem for a given starting state.
Two different versions with increasing levels of performance.
Each version has identical inputs and outputs.

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

module iteration_functions
export Solve_HH_Problem_v1, Solve_HH_Problem_v2
using LinearAlgebra

# -----------------------------------------------------------------------------------------------------
# Version 1 - using only base functions + for-loops.
# -----------------------------------------------------------------------------------------------------
function Solve_HH_Problem_v1(Value_Function, kt0_index, zt_index, params)

    # Unpack utility parameters and grids.
    α, β, δ = params.α, params.β, params.δ;
    k_values = params.k_values;
    z_values = params.z_values;
    z_probs = params.z_probs;

    # Get starting capital and productivity values from index.
    kt0 = k_values[kt0_index];
    zt = z_values[zt_index];

    # Variables to store candidate optimal values for the Value Function and Policy Function.
    v_max = -Inf;
    kt1_optimal = 0.0;

    # Check all possible next period capital choices.
    for kt1_index in eachindex(k_values)

        # Get capital value from index.
        kt1 = k_values[kt1_index];

        # Calculate the Value Function for given starting capital and next period capital choice.
        new_v_max = log(zt * (kt0 ^ α) + (1 - δ) * kt0 - kt1) + β * dot(Value_Function[kt1_index, :], z_probs[zt_index, :]);

        # Check if this capital choice gives the highest Value Function value.
        if new_v_max > v_max
        
            # Update candidate values.
            v_max = new_v_max;
            kt1_optimal = kt1;
        end
    end

    return v_max, kt1_optimal;
end

# -----------------------------------------------------------------------------------------------------
# Version 2 - using broadcast/vectorized calculation instead of for-loop.
# -----------------------------------------------------------------------------------------------------
function Solve_HH_Problem_v2(Value_Function, kt0_index, zt_index, params)

    # Unpack utility parameters and grids.
    α,β,δ = params.α, params.β, params.δ;
    k_values = params.k_values;
    z_values = params.z_values;
    z_probs = params.z_probs;

    # Get capital and productivity values from index.
    kt0 = k_values[kt0_index];
    zt = z_values[zt_index];

    # Calculate array of value function values for all next period capital choices.
    V_max_values = log.(zt * (kt0 ^ α) + (1 - δ) * kt0 .- k_values) + β * (Value_Function * z_probs[zt_index, :]);

    # Get index for the optimal capital choice.
    kt1_index_optimal = argmax(V_max_values);

    # Get the optimal Value Function and Policy Function values.
    kt1_optimal = k_values[kt1_index_optimal];
    v_max = V_max_values[kt1_index_optimal];

    return v_max, kt1_optimal;
end


# Benchmarking
#using BenchmarkTools
#@btime Solve_HH_Problem_v1(Value_Function, kt0_index, zt_index, params); # 49.3 μs
#@btime Solve_HH_Problem_v2(Value_Function, kt0_index, zt_index, params); #  6.1 μs

end # end module