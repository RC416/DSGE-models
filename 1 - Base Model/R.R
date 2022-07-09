# Base Model implementation in R

# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 1.000
number_of_iterations = 1000

# Calculate the steady-state level of capital.
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ^ (1/(alpha-1))  

# Create a grid of capital values around steady-state (+/- 50%).
number_of_k_values = 201
k_low_pct = 0.50
k_high_pct = 1.50
k_values = seq(k_low_pct*k_steady, k_high_pct*k_steady, length.out=number_of_k_values)

# Initialize Value Function and Policy Function (as arrays).
Value_Function = array(0, c(number_of_iterations, number_of_k_values))
Policy_Function = array(0, c(number_of_iterations, number_of_k_values))

# Perform value function iteration.
for (iteration in 2:number_of_iterations)
  {
    for (kt0_index in 1:number_of_k_values)       # for each level of starting capital...
    {
      # Variables to store candidate optimal values for Value Function and Policy Function.
      v_max = -Inf
      kt1_optimal = 0
     
         for (kt1_index in 1:number_of_k_values)
         {
           # Get capital values from index.
           kt0 = k_values[kt0_index]
           kt1 = k_values[kt1_index]
           
           # Calculate value of Value Function for given starting capital and next period capital choice.
           new_value_function_value = log((kt0^alpha) + (1-delta)*kt0 - kt1) + beta*Value_Function[iteration-1, kt1_index]
           
           # Check if this capital choice gives highest Value Function value.
           if (new_value_function_value > v_max)
           {
             # Update candidate values.
             v_max = new_value_function_value
             kt1_optimal = kt1
           }
         }
     
     # Update Value Function and Policy Function with optimal values.
     Value_Function[iteration, kt0_index] = v_max
     Policy_Function[iteration, kt0_index] = kt1_optimal
     }
}


# Plot various iterations of the Value Function.
plot(c(min(k_values),max(k_values)), c(min(Value_Function),max(Value_Function)), type="l", col="white", # hidden values to establish plot size
     main = "Value Function", xlab="k", ylab="V(k)")
lines(k_values, Value_Function[1,])
for (iteration in seq(number_of_iterations/10,number_of_iterations, number_of_iterations/10))
{
  lines(k_values, Value_Function[iteration,])
}
legend("right", legend = c(1,seq(number_of_iterations/10,number_of_iterations, number_of_iterations/10)))


# Plot final Policy Function.
plot(k_values, k_values, type="l", lty=2, col="black",
     main = "Policy Function", xlab = "k", ylab = "g(k)")
lines(k_values, Policy_Function[number_of_iterations,], col="blue")
legend("topleft", legend = c("g(k)", "45-degree line"), col=c("blue","black"), lty=c(1,2))
