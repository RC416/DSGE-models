% Base Model implementation in Matlab

% Assign parameter values.
alpha = 0.400;
beta  = 0.987;
delta = 1.000;
number_of_iterations = 1000;

% Calculate the steady-state level of capital.
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ^ (1/(alpha-1));

% Create a range of capital values around steady-state (+/- 50%).
number_of_k_values = 201;
k_low_pct = 0.50;
k_high_pct = 1.50;
k_values = linspace(k_low_pct*k_steady, k_high_pct*k_steady, number_of_k_values);

% Initialize Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values);

% Perform value function iteration.
for iteration = 2:(number_of_iterations)
    
    for ind_kt0 = 1:number_of_k_values      % for each level of starting captial...
        
        % Variables to store candidate optimal values for Value Function and Policy Function.
        v_max = -inf;
        kt1_optimal = 0;

        for ind_kt1 = 1:number_of_k_values  % ...check all next period capital choices

            % Calculate Value Function for given starting capital and next period capital choice.
            new_value_function_value = log((k_values(ind_kt0)^alpha) - k_values(ind_kt1) + (1-delta)*k_values(ind_kt0))...
                + beta*Value_Function(iteration-1, ind_kt1);
            
            % Check if this capital choice gives highest Value Function value.
            if new_value_function_value > v_max
                
                % Update candidate values.
                v_max = new_value_function_value;
                kt1_optimal = k_values(ind_kt1);
            end
        end

        % Update Value Function and Policy Function with optimal values.
        Value_Function(iteration, ind_kt0) = v_max;
        Policy_Function(iteration, ind_kt0) = kt1_optimal;
    end
end

% Plot various iterations of the Value Function.
figure(1)
hold on
for iteration = 1:number_of_iterations/10:number_of_iterations
    plot(k_values, Value_Function(iteration,:))
end
hold off
xlabel('k')
ylabel('V(k)')
title('Value Function')
legend(string(1:number_of_iterations/10:number_of_iterations))

% Plot final Policy Function.
figure(2)
hold on
plot(k_values,Policy_Function(number_of_iterations,:))
plot(k_values,k_values, '--', Color='k')
hold off
xlabel('k')
ylabel('g(k)')
title('Policy Function')
legend('g(k)','45^o Line', 'Location', 'northwest')