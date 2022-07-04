#=
Base Model implementation in Julia.
=#

# Assign parameter values.
α = 0.400;
β = 0.987;
δ = 1.000;
number_of_iterations = 1000;

# Calculate the steady-state level of capital.
k_steady = ((1-β*(1-δ))/(α*β)) ^ (1/(α-1))

# Create a range of capital values around steady-state (+/- 50%).
number_of_k_values = 201;
k_low_pct = 0.50;
k_high_pct = 1.50;
k_values = range(k_low_pct*k_steady, k_high_pct*k_steady, length=number_of_k_values);

# Initialize Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values);

# Perform value function iteration.
for iteration in 2:number_of_iterations

    for kt0_index in eachindex(k_values)      # for each level of starting capital...
        
        # Calculate Value Function for given starting capital and next period capital choice.
        v_max = -Inf;
        kt1_optimal = 0;

        for kt1_index in eachindex(k_values)  # ...check all next period capital choices 
            
            # Get capital values from index.
            kt0 = k_values[kt0_index];
            kt1 = k_values[kt1_index];

            # Calculate value of Value Function for given starting capital and next period capital choice.
            new_value_function_value = log(kt0^α + (1-δ)*kt0 - kt1) + β*Value_Function[iteration-1, kt1_index];

            # Check if this capital choice gives highest Value Function value.
            if new_value_function_value > v_max

                # Update candidate values.
                v_max = new_value_function_value;
                kt1_optimal = kt1;
            end
        end
        
        # Update Value Function and Policy function with optimal values.
        Value_Function[iteration, kt0_index] = v_max;
        Policy_Function[iteration, kt0_index] = kt1_optimal;
    end
end

# Plot various iterations of the Value Function.
using Plots
Figure1 = plot(k_values, Value_Function[1,:], label="1")
for iteration = number_of_iterations/10:number_of_iterations/10:number_of_iterations
    plot!(k_values, Value_Function[Int(iteration),:], label=string(Int(iteration)))
end
xlabel!("k")
ylabel!("V(k)")
title!("Value Function")
display(Figure1)

# Plot final Policy Function.
Figure2 = plot(k_values, Policy_Function[number_of_iterations,:], label="g(k)", legend=:topleft)
plot!(k_values,k_values, label="45-degree line", color=:black)
xlabel!("k")
ylabel!("g(k)")
title!("Policy Function")
display(Figure2)