% Function 1: Solve Household Problem (version 1)
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
    Solve_HH_Problem_v1(q, a_start_index, e_start_index, Value_Function, params)

% Unpack utility parameters and grids.
sigma = params.sigma;
beta = params.beta;
a_grid = params.a_grid;
e_grid = params.e_grid;
e_probs = params.e_probs;
number_of_a_values = params.number_of_a_values;

% Get value of state variables.
a_start = a_grid(a_start_index);
e_start = e_grid(e_start_index);

% Variables to store candidate optimal values for Value Function and Policy Function.
v_max = -Inf;
a_next_optimal = 0.0;
a_next_optimal_index = 0;

% Search over possible next period borrowing choices.
for a_next_index = 1:number_of_a_values
    
    % Get next credit value and the value of consumption implied by the budget constraint.
    a_next = a_grid(a_next_index);
    consumption = a_start + e_start - q*a_next;

    % Check budget constraint: if consumption is negative, skip this value.
    if (consumption <= 0); break; end
    
    % Calculate value function for given choice of next period capital.
    new_v_max = ((consumption)^(1-sigma))/(1-sigma) + beta*dot(Value_Function(a_next_index,:), e_probs(e_start_index,:));

    % Check if this capital choice gives highest Value Function value.
    if new_v_max > v_max

        % Update candidate values.
        v_max = new_v_max;
        a_next_optimal = a_next;
        a_next_optimal_index = a_next_index;
    end
end
end