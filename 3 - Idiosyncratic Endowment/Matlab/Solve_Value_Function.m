% Function 2: Solve Value Function
% 
% Solves for the Value Function and Policy Function using value function iteration.
% Applies the Solve Household Problem function to all possible starting states for each iteration.
% Search parameters (max iterations, tolerance, etc.) are defined in the function.
% 
% Input:
%     - Market price for bond
%     
% Output:
%     - Value Function
%     - Policy Function
%     - Policy Function (with indices instead of values)

function [Value_Function, Policy_Function, Policy_Function_Index] = ...
    Solve_Value_Function(q, params)

% Unpack relevant parameters.
number_of_a_values = params.number_of_a_values;
number_of_e_values = params.number_of_e_values;

% Arrays to hold 2 value function iterations.
Value_Function = zeros(number_of_a_values, number_of_e_values);
Policy_Function = zeros(number_of_a_values, number_of_e_values);
Policy_Function_Index = zeros(number_of_a_values, number_of_e_values);

Value_Function_New = zeros(number_of_a_values, number_of_e_values);
Policy_Function_New = zeros(number_of_a_values, number_of_e_values);
Policy_Function_Index_New = zeros(number_of_a_values, number_of_e_values);

% Iteration parameters.
dist = Inf;
iteration_count = 0;
max_iterations = 5000;
tolerance = 1e-6;

% Solve for the Value Function and Policy Function.
while (dist > tolerance) & (iteration_count < max_iterations)

    % Loop over all possible starting states.
    for a_start_index = 1:number_of_a_values
        for e_start_index = 1:number_of_e_values
                        
            % Solve the Value Function and Policy Function and update values.
            [V, g, g_index] = Solve_HH_Problem_v1(q, a_start_index, e_start_index, Value_Function, params);
            %[V, g, g_index] = Solve_HH_Problem_v2(q, a_start_index, e_start_index, Value_Function, params);

            Value_Function_New(a_start_index, e_start_index) = V;
            Policy_Function_New(a_start_index, e_start_index) = g;
            Policy_Function_Index_New(a_start_index, e_start_index) = g_index;
        end
    end
    
    % Update search parameters.
    dist = max(abs(Value_Function - Value_Function_New), [], 'all');
    iteration_count = iteration_count + 1;

    % Update the Value Function and Policy Function.
    Value_Function = Value_Function_New;
    Policy_Function = Policy_Function_New;
    Policy_Function_Index = Policy_Function_Index_New;
    
    % Print warning if convergence is not achieved.
    if iteration_count >= max_iterations
        fprintf("Warning: value function did not converge after %d iterations \n", iteration_count);
    end
end
end