# Function to perform an iteration of value function iteration.
# Two different versions with increasing levels of performance.
# Each version has identical inputs and outputs.
# 
# Input:
#   - Current value function
# - Starting capital and productivity level 
# 
# Output:
#   - Next value function value
# - Optimal choice of next period capital
# 
# Version
# v1: using only base functions + for-loops
# v2: using broadcast/vectorized calculation instead of for-loop

# -----------------------------------------------------------------------------------------------------
# Version 1 - using only base functions + for-loops.
# -----------------------------------------------------------------------------------------------------
Iterate_Value_Function_v1 = function(Previous_Value_Function, kt0_index, zt_index, params)
  {
  # Unpack utility parameters and grids.
  alpha = params$alpha
  beta = params$beta
  delta = params$delta
  k_values = params$k_values
  z_values = params$z_values
  z_probs = params$z_probs
  
  # Get starting capital and productivity values from index.
  kt0 = k_values[kt0_index]
  zt = z_values[zt_index]
  
  # Variables to store candidate optimal values for Value Function and Policy Function.
  v_max = -Inf
  kt1_optimal = 0
  
  for (kt1_index in 1:number_of_k_values)
  {
    # Get capital value from index.
    kt1 = k_values[kt1_index]
    
    # Calculate value function for given choice of next period capital.
    #new_value_function_value = log(zt*(kt0^alpha) + (1-delta)*kt0 - kt1) 
    #                            + beta*(Previous_Value_Function[kt1_index,]%*%z_probs[zt_index,])
    new_value_function_value = log(zt*(kt0^alpha) + (1-delta)*kt0 - kt1) + beta*sum(Previous_Value_Function[kt1_index,]*z_probs[zt_index,])
    
    # Check if this capital choice gives highest Value Function value.
    if (new_value_function_value > v_max)
    {
      # Update candidate values.
      v_max = new_value_function_value
      kt1_optimal = kt1
    }
  }
  return(list(v_max=v_max, kt1_optimal=kt1_optimal))
}
  
# -----------------------------------------------------------------------------------------------------
# Version 2 - using broadcast/vectorized calculation instead of for-loop.
# -----------------------------------------------------------------------------------------------------
Iterate_Value_Function_v2 = function(Previous_Value_Function, kt0_index, zt_index, params)
  {
  # Unpack utility parameters and grids.
  alpha    = params$alpha
  beta     = params$beta
  delta    = params$delta
  k_values = params$k_values
  z_values = params$z_values
  z_probs  = params$z_probs
  
  # Get capital and productivity values from index.
  kt0 = k_values[kt0_index]
  zt = z_values[zt_index]
  
  # Calculate array of value function values for all k_values.
  V_max_values = log(zt*(kt0^alpha) + (1-delta)*kt0 - k_values) + 
    beta*drop(Previous_Value_Function %*% z_probs[zt_index,])
  
  # Get the index for the optimal capital choice.
  kt1_index_optimal = which.max(V_max_values)
  
  # Get Value Function and Policy Function values.
  kt1_optimal = k_values[kt1_index_optimal]
  v_max = V_max_values[kt1_index_optimal]
  
  return(list(v_max=v_max, kt1_optimal=kt1_optimal))
}