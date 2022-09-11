% Stochastic Growth model implemented in Matlab.
% 
% Steps:
%     1 - Define utility parameters, grids, and parameter struct.
%     2 - Perform value function iteration.
%     3 - Plot results.

%% 1 - Define utility parameters, grids, and parameter struct.

% Assign parameter values.
alpha = 0.400;
beta  = 0.987;
delta = 0.012;
number_of_iterations = 1000;

% Calculate the steady-state level of capital.
k_steady = ((1 - beta * (1 - delta)) / (alpha * beta)) ^ (1 / (alpha - 1));

% Create a range of capital values around the steady-state (+/- 2%).
number_of_k_values = 201;
k_low_pct = 0.98;
k_high_pct = 1.02;
k_values = linspace(k_low_pct * k_steady, k_high_pct * k_steady, number_of_k_values)';

% Get productivity levels and transition probabilities.
z_probs = readmatrix("Inputs\z_probs.csv");
z_values = readmatrix("Inputs\z_values.csv");
number_of_z_values = size(z_values, 1);

% Initialize the Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_k_values, number_of_z_values, number_of_iterations);
Policy_Function = zeros(number_of_k_values, number_of_z_values, number_of_iterations);

% Store utility parameters and capital/productivity grids in a struct for passing to a function.
params = struct( ...
    alpha=alpha, beta=beta, delta=delta, ...
    k_values=k_values, z_values=z_values, ...
    z_probs=z_probs);

%% 2 - Perform value function iteration.

for iteration = 2:number_of_iterations

    % Loop over all possible starting states.
    for kt0_index = 1:number_of_k_values
        for zt_index = 1:number_of_z_values

            % Solve the Value Function and Policy Function and update values.
            %[V, g] = Solve_HH_Problem_v1(Value_Function(:, :, iteration - 1), kt0_index, zt_index, params);
            [V, g] = Solve_HH_Problem_v2(Value_Function(:, :, iteration - 1), kt0_index, zt_index, params);
    
            Value_Function(kt0_index, zt_index, iteration) = V;
            Policy_Function( kt0_index, zt_index, iteration) = g;
        end
    end
end

%% 3 - Plot results.

% Plot the Value Function for different starting states.
figure(1)
hold on
for zt_index = 1:number_of_z_values
    plot(k_values,Value_Function(:, zt_index, number_of_iterations))
end
hold off
xlabel('k')
ylabel('V(k,z)')
title('Value Function')
lgd = legend([string(round(flip(z_values), 2))]);
lgd.Title.String = "z values";

% Plot the final Policy Function for certain productivity values.
figure(2)
z_indices = [1, 4, 6, 8, 11];
hold on
for z_index = 1:size(z_indices, 2)
    plot(k_values,Policy_Function(:, z_indices(z_index), number_of_iterations)')
end
plot(k_values, k_values, '--', Color='k')
hold off
xlabel('k')
ylabel('g(k,z)')
title('Policy Function')
lgd = legend([string(round(z_values(z_indices), 2))], 'Location', 'northwest');
lgd.Title.String = "z values";