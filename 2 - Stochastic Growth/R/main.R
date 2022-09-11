# Stochastic Growth model implemented in R.
# 
# Steps:
#     1 - Define utility parameters, grids, and parameter struct.
#     2 - Create a function to solve the household's problem for a given starting state.
#     3 - Perform value function iteration.
#     4 - Plot results.

# -----------------------------------------------------------------------------------------------------
# 1 - Define utility parameters, grids, and parameter struct.
# -----------------------------------------------------------------------------------------------------
# Assign parameter values.
alpha = 0.400
beta  = 0.987
delta = 0.012
number_of_iterations = 1000

# Calculate the steady-state level of capital.
k_steady = ((1 - beta * (1 - delta)) / (alpha * beta)) ^ (1 / (alpha - 1)) 

# Create a range of capital values around steady-state (+/- 2%).
number_of_k_values = 201
k_low_pct = 0.98
k_high_pct = 1.02
k_values = seq(k_low_pct * k_steady, k_high_pct * k_steady, length.out=number_of_k_values)

# Get productivity levels and transition probabilities.
z_probs = as.matrix(read.csv("Inputs/z_probs.csv", sep=",", header=FALSE))
z_values = as.matrix(read.csv("Inputs/z_values.csv", sep=",", header=FALSE))
number_of_z_values = dim(z_values)[1]

# Initialize the Value Function and Policy Function (as arrays).
Value_Function = array(0, c(number_of_iterations, number_of_k_values, number_of_z_values))
Policy_Function = array(0, c(number_of_iterations, number_of_k_values, number_of_z_values))

# Store utility parameters and capital/productivity grids in a list for passing to a function.
params = list(alpha=alpha, beta=beta, delta=delta,
              k_values=k_values, z_values=z_values, z_probs=z_probs)

# -----------------------------------------------------------------------------------------------------
# 2 - Create a function to solve Value Function for given starting states.
# -----------------------------------------------------------------------------------------------------
Solve_HH_Problem = function(Value_Function, kt0_index, zt_index, params){
  
  # Unpack utility parameters and grids.
  alpha = params$alpha
  beta = params$beta
  delta = params$delta
  k_values = params$k_values
  z_values = params$z_values
  z_probs = params$z_probs
  
  # Get capital and productivity values from index.
  kt0 = k_values[kt0_index]
  zt = z_values[zt_index]
  
  # Calculate array of value function values for all capital values.
  V_max_values = log(zt * (kt0 ^ alpha) + (1 - delta) * kt0 - k_values) + 
                  beta * drop(Value_Function %*% z_probs[zt_index, ])
  

  # Get the index for the optimal capital choice.
  kt1_index_optimal = which.max(V_max_values)
  
  # Get the Value Function and Policy Function values.
  kt1_optimal = k_values[kt1_index_optimal]
  v_max = V_max_values[kt1_index_optimal]
  
  return(list(v_max=v_max, kt1_optimal=kt1_optimal))
}

# Alternative: import other versions from script.
# source("custom_functions")
# use Solve_HH_Problem_v1 or Solve_HH_Problem_v2

# -----------------------------------------------------------------------------------------------------
# 3 - Perform value function iteration.
# -----------------------------------------------------------------------------------------------------
for (iteration in 2:number_of_iterations)
{
  # Loop over all possible starting states.
  for (kt0_index in 1:number_of_k_values)
  {
    for (zt_index in 1:number_of_z_values)
    {
      # Solve the Value Function and Policy Function and update values.
      result = Solve_HH_Problem(Value_Function[iteration-1, , ], kt0_index, zt_index, params)
      
      Value_Function[iteration,kt0_index,zt_index] = result$v_max
      Policy_Function[iteration,kt0_index,zt_index] = result$kt1_optimal
    }
  }
}

# -----------------------------------------------------------------------------------------------------
# 4 - Plot results.
# ----------------------------------------------------------------------------------------------------- 
# Plot the Value Function for different starting states.
plot(c(min(k_values), max(k_values)), 
     c(min(Value_Function[number_of_iterations, , ]), (Value_Function[number_of_iterations, , ])),
     type="l", col="white", # hidden values to establish plot size
     main="Value Function", xlab="k", ylab="V(k,z)")
for (z_index in 1:number_of_z_values)
{
  lines(k_values, Value_Function[number_of_iterations, ,z_index])
}
  legend("right", legend=round(rev(z_values), 2), title="z values")

# Plot the final Policy Function for certain productivity values.
z_indices = c(1, 4, 6, 8, 11)
plot(k_values, k_values, type="l", lty=2, col="black",
     main="Policy Function", xlab="k", ylab="g(k)")
for (z_index in z_indices){
  lines(k_values, Policy_Function[number_of_iterations, , z_index])
}
legend("topleft", legend=c(round(rev(z_values[z_indices]), 2), "45-degree line"),
       lty=c(1, 1, 1, 1, 1, 2), title="z values")