% Function 3: Get Population Distribution. 
% 
% Solves for the steady state distribution over credit and endowment states given
% a policy function (with index values).
% Search parameters (max iterations, tolerance, etc.) are defined in the function.
% 
% Input:
%     - Policy Function (with index values)
% 
% Output:
%     - Steady-state population distribution

function Population_Distribution = ...
    Get_Population_Distribution(Policy_Function_Index, params)

% Unpack relevant parameters.
e_probs = params.e_probs;
number_of_a_values = params.number_of_a_values;
number_of_e_values = params.number_of_e_values;

% Arrays to store 2 iterations of finding population distribution.    
Population_Distribution = ones(number_of_a_values, number_of_e_values) ./ (number_of_a_values * number_of_e_values);
New_Distribution = Population_Distribution;

% Iteration parameters.
dist = Inf;
iteration_count = 0;
max_iterations = 5000;
tolerance = 1e-10;

% Solve for population distribution.
while (dist > tolerance) & (iteration_count < max_iterations)

    % Get "inflow" to each credit-endowment state in next period.
    for a_index = 1:number_of_a_values
        for e_index = 1:number_of_e_values
            
            % Sum distribution-weighted inflow into given state.
            inflow = sum( (Population_Distribution .* (Policy_Function_Index == a_index)) * e_probs(:,e_index) );
            New_Distribution(a_index, e_index) = inflow;
        end
    end

    % Update search parameters.
    dist = sum(abs(Population_Distribution - New_Distribution), 'all');
    iteration_count = iteration_count + 1;

    % Update population distribution.
    Population_Distribution = New_Distribution;
    
    % Print warning if convergence is not achieved.
    if iteration_count >= max_iterations
        fprintf("Warning: population distribution did not converge after %d iterations \n", iteration_count);
    end
end
end