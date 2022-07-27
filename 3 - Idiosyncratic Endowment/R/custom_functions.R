# Three functions with multiple versions:
#   1. Solve Household Problem
#      - v1, v2
#   2. Solve Value Function
#      - Uses Solve Household Problem to solve for the value function and policy function
#   3. Get Population Distribution

# -----------------------------------------------------------------------------------------------------
# Function 1 - Version 1 - Solve Household Problem using only base functions + for-loops.
# -----------------------------------------------------------------------------------------------------

# Function find the next iteration of the Value Function and Policy Function 
# by solving the household's problem for a given starting state.
# Four different versions with increasing levels of performance.
# Each version has identical inputs and outputs.
# Input:
#    - Current Value Function
#    - Market price for bond
#    - Starting wealth and endowment level 

# Output:
#    - Next value function value
#    - Optimal choice of savings/borrowing for next period

# Version
#    v1: using only for-loops
#    v2: using arrays / broadcast instead of for-loop

Solve_HH_Problem_v1 = function(q, a_start_index, e_start_index, Value_Function, params)
{
  # Unpack utility parameters and grids.
  sigma   = params$sigma
  beta    = params$beta
  a_grid  = params$a_grid
  e_grid  = params$e_grid
  e_probs = params$e_probs
  
  # Get value of state variables.
  a_start = a_grid[a_start_index]
  e_start = e_grid[e_start_index]
  
  # Variables to store candidate optimal values for Value Function and Policy Function.
  v_max = -Inf
  a_next_optimal = 0.0
  a_next_optimal_index = 0
  
  # Search over possible next period borrowing choices.
  for (a_next_index in 1:number_of_a_values)
  {
    # Get next credit value and the value of consumption implied by the budget constraint.
    a_next = a_grid[a_next_index]
    consumption = a_start + e_start - q*a_next
    
    # Check budget constraint: if consumption is negative, skip this value.
    if (consumption <= 0) { next }
    
    # Calculate value function for given choice of next period capital.
    new_v_max = ((consumption)^(1-sigma))/(1-sigma) + beta*sum(Value_Function[a_next_index, ] * e_probs[e_start_index, ]);
    
    # Check if this capital choice gives highest Value Function value.
    if (new_v_max > v_max)
    {
      # Update candidate values.
      v_max = new_v_max
      a_next_optimal = a_next
      a_next_optimal_index   = a_next_index
    }
  }
  return(list(v_max = v_max,
              a_next_optimal = a_next_optimal,
              a_next_optimal_index = a_next_optimal_index))
}

# -----------------------------------------------------------------------------------------------------
# Function 1 - Version 2 - Solve Household Problem using vectorized calculation.
# -----------------------------------------------------------------------------------------------------
Solve_HH_Problem_v2 = function(q, a_start_index, e_start_index, Value_Function, params)
{
  # Unpack utility parameters and grids.
  sigma   = params$sigma
  beta    = params$beta
  a_grid  = params$a_grid
  e_grid  = params$e_grid
  e_probs = params$e_probs
  
  # Get value of state variables.
  a_start = a_grid[a_start_index]
  e_start = e_grid[e_start_index]
  
  # Vector of consumption values dictated by credit selection.
  Consumption = a_start + e_start - q*a_grid
  valid_indices = (Consumption > 0)
  
  # Calculate value function values.
  V_max_values = Consumption[valid_indices]^(1-sigma) / (1-sigma) + beta*(Value_Function[valid_indices, ] %*% e_probs[e_start_index, ])
  
  # Get index of optimal credit choice within the vector of valid indices.
  optimal_subindex = which.max(V_max_values)
  
  # Get values and original index of optimal value.
  a_next_optimal = (a_grid[valid_indices])[optimal_subindex]
  v_max = V_max_values[optimal_subindex]
  a_next_optimal_index = match(a_next_optimal, a_grid)
  
  return(list(v_max = v_max,
              a_next_optimal = a_next_optimal,
              a_next_optimal_index = a_next_optimal_index))
}

# -----------------------------------------------------------------------------------------------------
# Function 2 - Solve for Value Function and Policy Function.
# -----------------------------------------------------------------------------------------------------
# Solves for the Value Function and Policy Function using value function iteration.
# Applies Solve Household Problem to all possible starting states for each iteration.
# Search parameters (max iterations, tolerance, etc.) are defined in the function.

# Input:
#  - Market price for bond

# Output:
#  - Value Function
#  - Policy Function
#  - Policy Function (with indices instead of values)

Solve_Value_Function = function(q, params)
{
  # Unpack relevant parameters.
  a_grid  = params$a_grid
  e_grid  = params$e_grid
  number_of_a_values = params$number_of_a_values
  number_of_e_values = params$number_of_e_values
  
  # Arrays to hold 2 value function iterations.
  Value_Function = array(0, c(number_of_a_values, number_of_e_values))
  Policy_Function = array(0, c(number_of_a_values, number_of_e_values))
  Policy_Function_Index = array(0, c(number_of_a_values, number_of_e_values))
  
  Value_Function_New = array(0, c(number_of_a_values, number_of_e_values))
  Policy_Function_New = array(0, c(number_of_a_values, number_of_e_values))
  Policy_Function_Index_New = array(0, c(number_of_a_values, number_of_e_values))
  
  # Iteration parameters.
  dist = Inf
  iteration_count = 0
  max_iterations = 5000
  tolerance = 1e-6
  
  # Solve for Value Function and Policy Function.
  while ((dist > tolerance) & (iteration_count < max_iterations))
  {
    # Loop over all possible starting states.
    for (a_start_index in 1:number_of_a_values)
    {
      for (e_start_index in 1:number_of_e_values)
      {
        # Solve Value Function and Policy Function and update values.
        #result = Solve_HH_Problem_v1(q, a_start_index, e_start_index, Value_Function, params);
        result  = Solve_HH_Problem_v2(q, a_start_index, e_start_index, Value_Function, params);
        
        Value_Function_New[a_start_index, e_start_index] = result$v_max;
        Policy_Function_New[a_start_index, e_start_index] = result$a_next_optimal;
        Policy_Function_Index_New[a_start_index, e_start_index] = result$a_next_optimal_index;
      }
    }
    
    # Update search parameters.
    dist = max(abs(Value_Function - Value_Function_New))
    iteration_count = iteration_count + 1
    
    # Update Value Function.
    Value_Function = Value_Function_New
    Policy_Function = Policy_Function_New
    Policy_Function_Index = Policy_Function_Index_New
    
    # Print warning if convergence is not achieved.
    if (iteration_count >= max_iterations)
    {
      print(sprintf("Warning: value function did not converge after %s iterationts", iteration_count))
    }
  }
  return(list(Value_Function = Value_Function,
              Policy_Function = Policy_Function,
              Policy_Function_Index = Policy_Function_Index))
}

# -----------------------------------------------------------------------------------------------------
# Function 3 - Get Population Distribution from Policy Function.
# -----------------------------------------------------------------------------------------------------
# Solves for the steady state distribution over credit and endowment states given
# a policy function (with index values).
# Search parameters (max iterations, tolerance, etc.) are defined in the function.

# Input:
#   - Policy Function (with index values)

# Output:
#   - Steady-state population distribution

Get_Population_Distribution = function(Policy_Function_Index, params)
{
  # Unpack relevant parameters.
  a_grid = params$a_grid
  e_grid = params$e_grid
  e_probs = params$e_probs
  number_of_a_values = params$number_of_a_values
  number_of_e_values = params$number_of_e_values
  
  # Arrays to store 2 iterations of finding population distribution.    
  Population_Distribution = array(1, c(number_of_a_values, number_of_e_values)) / (number_of_a_values * number_of_e_values)
  New_Distribution = Population_Distribution
  
  # Iteration parameters.
  dist = Inf
  iteration_count = 0
  max_iterations = 5000
  tolerance = 1e-10
  
  # Solve for population distribution.
  while ((dist > tolerance) & (iteration_count < max_iterations))
  {
    # Get "inflow" to each credit-endowment state in next period.
    for (a_index in 1:number_of_a_values)
    {
      for (e_index in 1:number_of_e_values)
      {
        # Sum distribution-weighted inflow into given state.
        inflow = sum( (Population_Distribution * (Policy_Function_Index == a_index)) %*% e_probs[ ,e_index] )
        New_Distribution[a_index, e_index] = inflow
      }
    }
    
    # Update search parameters.
    dist = sum(abs(Population_Distribution - New_Distribution))
    iteration_count = iteration_count + 1
    
    # Update population distribution.
    Population_Distribution = New_Distribution
    
    # Print warning if convergence is not achieved.
    if (iteration_count >= max_iterations)
    {
      print(sprintf("Warning: population distribution did not converge after %s iterations", iteration_count))
    }
  }
  return(Population_Distribution)
}