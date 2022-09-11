#=
Idiosyncratic Endowment model implemented in Julia.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Solve model and get market bond price using binary search.
        - Guess bond price
        - Solve for the Value Function and Policy Function
        - Get the distribution of credit and productivity levels
        - Check the market clearing condition and update bond price guess
    3 - Plot results.
=#

include("custom_functions.jl")
using .custom_functions

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------

# Endowment parameters.
e_high = 1.0;                                               # high endowment
e_low  = 0.1;                                               # low endowment 
e_grid = [e_low, e_high];                                   # endowment grid
number_of_e_values = 2;                                     
e_probs = [0.500 0.500 ; 0.075 0.925];                      # transition probabilities

# Utility parameters.
σ = 1.5;                                                    # risk aversion coefficient
β = 0.99322;                                                # discount factor

# Credit parameters.
a_high = 4;                                                 # upper credit limit
a_low  = -2;                                                # lower credit limit / borrowing constraint
number_of_a_values = 100;                                   # credit grid size
a_grid = LinRange(a_low, a_high, number_of_a_values);       # credit grid

# Store parameters in a struct for passing to a function.
struct Parameters
    σ::Float64
    β::Float64
    a_grid::Array{Float64, 1}
    e_grid::Array{Float64, 1}
    e_probs::Array{Float64, 2}
    number_of_a_values::Int64
    number_of_e_values::Int64
end

params = Parameters(σ, β, a_grid, e_grid, e_probs, number_of_a_values, number_of_e_values);

# -----------------------------------------------------------------------------------------------------
# 2 - Solve model and get market bond price using binary search.
# -----------------------------------------------------------------------------------------------------

# Range of bond price values to search.
q_min = 0.985;
q_max = 1.100;

# Optional: floor for q_min such that those at credit limit can still afford positive consumption.
q_min = (a_low + e_low) / a_low;

# Placeholder for the market clearing condition.
mcc = Inf;

# Iteration parameters.
dist = Inf;
iteration_count = 0;
max_iterations = 20;
tolerance = 1e-3;

# Solve for market price q.
while (dist > tolerance) & (iteration_count < max_iterations)

    # Get value of q from middle of range.
    q = (q_min + q_max)/2;

    # Solve for the Value Function and Policy Function.
    global Value_Function, Policy_Function, Policy_Function_Index = Solve_Value_Function(q, params);

    # Get the Population Distribution.
    global Population_Distribution = Get_Population_Distribution(Policy_Function_Index, params);

    # Check the market clearing condition.
    global mcc = sum(Policy_Function .* Population_Distribution);

    # Update search parameters.
    global dist = abs(mcc);
    global iteration_count += 1;

    # Update range of q according to the sign of the market clearing condition.
    if mcc > 0; global q_min = q; end
    if mcc < 0; global q_max = q; end

    # Print results.
    println("Iteration $iteration_count: q=$(round(q, digits=6)), mcc=$(round(mcc,digits=6))");
    if iteration_count >= max_iterations; println("Warning: search for q did not converge after $iteration_count iterations"); end

end

# -----------------------------------------------------------------------------------------------------
# 3 - Plot results.
# ----------------------------------------------------------------------------------------------------- 
using Plots

# Policy functions (Figure 1. from Huggett 1993).
Figure1 = plot(legend=:bottomright)
plot!(a_grid, Policy_Function[:, 2], label="high entitlement")
plot!(a_grid, Policy_Function[:, 1], label="low entitlement")
plot!(a_grid, a_grid, linestyle=:dash, label="45-degree line", color=:black)
xlabel!("starting credit level")
ylabel!("optimal new credit level")
title!("Policy Function")
display(Figure1)

# Distribution of credit levels (Figure 2. from Huggett 1993).
Figure2 = plot(xlim=[-2, 1], ylim=[0, 1], legend=:topleft)
plot!(a_grid, cumsum(Population_Distribution[:, 2]), label="high entitlement")
plot!(a_grid, cumsum(Population_Distribution[:, 1]), label="low entitlement")
xlabel!("starting credit level")
title!("Cumulative Distribution Function for Credit Level")