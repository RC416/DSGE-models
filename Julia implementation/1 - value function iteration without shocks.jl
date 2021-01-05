

# assign initial values
α = 0.400
β = 0.987
δ = 0.012
ρ = 0.950
σ = 0.007
iterations = 100

zt = 1

# calculate steady-state level of capital
k_steady = ((1-β*(1-δ))/(α*β)) ^ (1/(α-1))                 # 100.44

# create deviations around steady state
k_values = range(0.98*k_steady, 1.02*k_steady, length=21)  # -2% to +2% of k_steady
k_values = round.(k_values, digits=2)                      # broadcast .



