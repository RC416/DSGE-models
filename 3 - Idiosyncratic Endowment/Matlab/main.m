% Idiosyncratic Entitlement model implemented in Matlab.
% 
% Steps:
%     1 - Define utility parameters, grids, and parameter struct.
%     2 - Solve model and get market bond price using binary search.
%         - Guess bond price
%         - Solve for Value Function and Policy Function
%         - Get distribution of credit and productivity levels
%         - Check market clearing condition and update bond price guess
%     3 - Plot results.

%% 1 - Define utility parameters, grids, and parameter struct.
tic
% Endowment parameters.
e_high = 1.0;                                               % high entitlement
e_low  = 0.1;                                               % low entitlement 
e_grid = [e_low, e_high];                                   % entitlement grid
number_of_e_values = 2;                                     
e_probs = [0.500 0.500 ; 0.075 0.925];                      % transition probabilities

% Utility parameters.
sigma = 1.5;                                                % risk aversion coefficient
beta = 0.99322;                                             % discount factor

% Credit parameters.
a_high = 4;                                                 % upper credit limit
a_low  = -2;                                                % lower credit limit / borrowing constraint
number_of_a_values = 500;                                   % credit grid size
a_grid = linspace(a_low, a_high, number_of_a_values);       % credit grid

% Store parameters in a struct for passing to a function.
params = struct(...
    sigma=sigma, beta=beta, ...
    a_grid=a_grid, e_grid=e_grid, e_probs=e_probs, ...
    number_of_a_values=number_of_a_values, ...
    number_of_e_values=number_of_e_values);

%% 2 - Solve model and get market bond price using binary search.

% Range of bond price values to search.
q_min = 0.985;
q_max = 1.100;

% Optional: floor for q_min such that those at credit limit can still afford positive consumption.
q_min = (a_low + e_low) / a_low;

% Placeholder for market clearing condition.
mcc = Inf;

% Iteration parameters.
dist = Inf;
iteration_count = 0;
max_iterations = 20;
tolerance = 1e-3;

% Solve for market price q.
while (dist > tolerance) & (iteration_count < max_iterations)

    % Get value of q from middle of range.
    q = (q_min + q_max)/2;

    % Solve for Value Function and Policy Function.
    [Value_Function, Policy_Function, Policy_Function_Index] = Solve_Value_Function(q, params);

    % Get Population Distribution.
    Population_Distribution = Get_Population_Distribution(Policy_Function_Index, params);

    % Check market clearing condition.
    mcc = sum(Policy_Function .* Population_Distribution, 'all');

    % Update search parameters.
    dist = abs(mcc);
    iteration_count = iteration_count + 1;

    % Update range of q according to the sign of mcc.
    if mcc > 0; q_min = q; end
    if mcc < 0; q_max = q; end

    % Print results.
    fprintf("Iteration %d: q=%f, mcc=%f \n", iteration_count, round(q,6), round(mcc,6));
    if iteration_count >= max_iterations; fprintf("Warning: search for q did not converge after %d iterations \n", iteration_count); end
end

%% Plot Results

% Policy functions (Figure 1. from Hugget 1993).
figure(1)
hold on
plot(a_grid, Policy_Function(:,2))
plot(a_grid, Policy_Function(:,1))
plot(a_grid, a_grid, '--', Color='k')
hold off
xlabel("starting credit level")
ylabel("optimal new credit level")
title("Policy Function")
legend(["high endowment", "low endowment", "45-degree line"], 'Location', 'southeast')

% Distribution of credit levels (Figure 2. from Huggett 1993).
figure(2)
hold on
plot(a_grid, cumsum(Population_Distribution(:,2)))
plot(a_grid, cumsum(Population_Distribution(:,1)))
hold off
xlim([-2,1])
ylim([0,1])
xlabel("starting credit level")
title("Cumulative Distribution Function for Credit Level")
legend(["high endowment", "low endowment"], 'Location', 'northwest')

toc