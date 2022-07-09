% Function to perform an iteration of value function iteration.
% Two different versions with increasing levels of performance.
% Each version has identical inputs and outputs.
% This version uses Matlab's vectorized log function.
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

function [v_max, kt1_optimal] = Iterate_Value_Function_v2(Previous_Value_Function, kt0_index, zt_index, params)

% Unpack utility parameters and grids.
alpha = params.alpha;
beta = params.beta;
delta = params.delta;
k_values = params.k_values;
z_values = params.z_values;
z_probs = params.z_probs;

% Get capital and productivity values from index.
kt0 = k_values(kt0_index);
zt = z_values(zt_index);

% Calculate array of value function values for all k_values.
V_max_values = log(zt*(kt0^alpha) + (1-delta)*kt0 - k_values) + beta*(Previous_Value_Function*z_probs(zt_index,:)');

% Get Value Function and Policy Function values.
[v_max, kt1_index_optimal] = max(V_max_values);
kt1_optimal = k_values(kt1_index_optimal);
end

