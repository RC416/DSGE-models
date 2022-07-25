# Idiosyncratic Entitlement model implemented in R.

# Steps:
#   1 - Define utility parameters, grids, and parameter struct.
#   2 - Solve model and get market bond price using binary search.
#       - Guess bond price
#       - Solve for Value Function and Policy Function
#       - Get distribution of credit and productivity levels
#       - Check market clearing condition and update bond price guess
#   3 - Plot results.

source("custom_functions.R")

# Endowment parameters.
e_high = 1.0                                                # high entitlement
e_low  = 0.1                                                # low entitlement 
e_grid = c(e_low, e_high)                                   # entitlement grid
number_of_e_values = 2                                      
e_probs = matrix(c(0.500, 0.500, 0.075, 0.925),
                 byrow=T, ncol=2)                           # transition probabilities

# Utility parameters.
sigma = 1.5                                                 # risk aversion coefficient
beta  = 0.99322                                             # discount factor

# Credit parameters.
a_high = 4                                                  # upper credit limit
a_low  = -2                                                 # lower credit limit / borrowing constraint
number_of_a_values = 50                                     # credit grid size
a_grid = seq(a_low, a_high, length.out=number_of_a_values)  # credit grid

# Store parameters in a named list for passing to a function.
params = list(sigma   = sigma,
              beta    = beta,
              a_grid  = a_grid,
              e_grid  = e_grid,
              e_probs = e_probs,
              number_of_a_values = number_of_a_values,
              number_of_e_values = number_of_e_values)

# -----------------------------------------------------------------------------------------------------
# 2 - Solve model and get market bond price using binary search.
# -----------------------------------------------------------------------------------------------------

# Range of bond price values to search.
q_min = 0.985
q_max = 1.100

# Optional: floor for q_min such that those at credit limit can still afford positive consumption.
q_min = (a_low + e_low) / a_low

# Placeholder for market clearing condition.
mcc = Inf;

# Iteration parameters.
dist = Inf
iteration_count = 0
max_iterations = 20
tolerance = 1e-3

# Solve for market price q.
while ((dist > tolerance) & (iteration_count < max_iterations))
{
  # Get value of q from middle of range.
  q = (q_min + q_max)/2
  
  # Solve for Value Function and Policy Function.
  vf_result = Solve_Value_Function(q, params)
  
  # Get Population Distribution.
  Population_Distribution = Get_Population_Distribution(vf_result$Policy_Function_Index, params)
  
  # Check market clearing condition.
  mcc = sum(vf_result$Policy_Function * Population_Distribution)
  
  # Update search parameters.
  dist = abs(mcc)
  iteration_count = iteration_count + 1
  
  # Update range of q according to the sign of mcc.
  if (mcc > 0) { q_min = q }
  if (mcc < 0) { q_max = q }
  
  # Print results.
  print(sprintf("Iteration %s: q=%f, mcc=%f", iteration_count, round(q, digits=6), round(mcc, digits=6)))
  if (iteration_count >= max_iterations)
    { 
    print(sprintf("Warning: search for q did not converge after %s iterations", iteration_count))
    }
}

# -----------------------------------------------------------------------------------------------------
# 3 - Plot results.
# ----------------------------------------------------------------------------------------------------- 

# Policy functions (Figure 1. from Huggett 1993).
plot(c(-2,1),c(-2,1), type="l", col="white",
     main="Policy Function", xlab="starting credit level", ylab="optimal new credit level")
lines(a_grid, vf_result$Policy_Function[ ,2], col="blue")
lines(a_grid, vf_result$Policy_Function[ ,1], col="green")
lines(a_grid, a_grid, lty=2, col="black")
legend("bottomright", legend = c("high endowment", "low endowment", "45-degree line"))

# Distribution of credit levels (Figure 2. from Huggett 1993).
plot(c(-2,1),c(0,1), type="l", col="white",
     main="Cumulative Distribution Function for Credit Level", xlab="starting credit level", ylab="")
lines(a_grid, cumsum(Population_Distribution[ ,2]), col="blue")
lines(a_grid, cumsum(Population_Distribution[ ,1]), col="green")
legend("topleft", legend = c("high endowment", "low endowment"), col=c("blue","green"))
