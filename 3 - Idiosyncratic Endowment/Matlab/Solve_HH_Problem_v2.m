% Function 1: Solve Household Problem (version 2)
%     
% Function to find the next iteration of the Value Function and Policy Function  
% by solving the household's problem for a given starting state.
% Two different versions with increasing levels of performance.
% Each version has identical inputs and outputs.
% 
% Input:
%     - Current Value Function
%     - Market price for bond
%     - Starting wealth and endowment level 
%     
% Output:
%     - Next value function value
%     - Optimal choice of savings/borrowing for next period
%     
% Version
%     v1: using only for-loops
%     v2: using arrays / broadcast instead of for-loop

function [v_max, a_next_optimal, a_next_optimal_index] = ...
    Solve_HH_Problem_v2(q, a_start_index, e_start_index, Value_Function, params)

% Unpack utility parameters and grids.
sigma = params.sigma;
beta = params.beta;
a_grid = params.a_grid;
e_grid = params.e_grid;
e_probs = params.e_probs;

% Get value of state variables.
a_start = a_grid(a_start_index);
e_start = e_grid(e_start_index);

% Vector of consumption values dictated by possible next period borrowing choices.
Consumption = a_start + e_start - q*a_grid;
valid_indices = (Consumption > 0);

% Calculate the Value Function values.
V_max_values = (Consumption(valid_indices).^(1-sigma))./(1-sigma) + beta*(Value_Function(valid_indices,:)*e_probs(e_start_index,:)')';

% Get the value and index of optimal value.
[v_max, optimal_subindex] = max(V_max_values);

% Get the values and original index of optimal value.
a_grid_valid = a_grid(valid_indices);
a_next_optimal = a_grid_valid(optimal_subindex);
a_next_optimal_index = find(a_grid == a_next_optimal);
end