% Function to solve the household's problem for a given starting state.
% Two different versions with increasing levels of performance.
% Each version has identical inputs and outputs.
% This version uses only for-loops.
% 
% Input:
%     - Current value function
%     - Starting capital and productivity level 
%     
% Output:
%     - Next value function value
%     - Optimal choice of next period capital
%     
% Version
%     v1: using only base functions + for-loops
%     v2: using broadcast/vectorized calculation instead of for-loop

function [v_max, kt1_optimal] = Solve_HH_Problem_v1(Value_Function, kt0_index, zt_index, params)

% Unpack utility parameters and grids.
alpha = params.alpha;
beta = params.beta;
delta = params.delta;
k_values = params.k_values;
z_values = params.z_values;
z_probs = params.z_probs;

% Get starting capital and productivity values from index.
kt0 = k_values(kt0_index);
zt = z_values(zt_index);

% Variables to store candidate optimal values for the Value Function and Policy Function.
v_max = -Inf;
kt1_optimal = 0.0;

% Check all possible next period capital choices.
for kt1_index = 1:size(k_values,1)
    
    % Get capital value from index.
    kt1 = k_values(kt1_index);

    % Calculate the Value Function for given starting capital and next period capital choice.
    new_v_max = log(zt * (kt0 ^ alpha) + (1 - delta) * kt0 - kt1) + ...
        beta * dot(Value_Function(kt1_index, :), z_probs(zt_index, :));

    % Check if this capital choice gives the highest Value Function value.
    if new_v_max > v_max
    
        % Update candidate values.
        v_max = new_v_max;
        kt1_optimal = kt1;
    end
end