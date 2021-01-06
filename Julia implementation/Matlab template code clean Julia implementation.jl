#=
Template code from Matlab implemented in Julia
=#

iter_max=160 # Choose number of iterations to make sure that V converges

α=0.3
β=0.65
k_bar=(α*β)^(1/(1-α))
# Make sure that your grid points for k include the steady state value of k
# K=0.05:0.025:0.15 
K=0.05:0.01:0.15 # finer grid points

# number of k values
N = length(K)

# initialize value function 
V = zeros(iter_max, N)
V[1,:]=zeros(1,N)  # initial guess 

# initialize capital policy function and temporary value function for the loop (W)
g = zeros(iter_max, N)
W = zeros(iter_max, N, N)

for t=2:iter_max
    for i=1:N
        vmax=-100000000
        for j=1:N
            W[t,i,j]=log(K[i]^α-K[j])+β*V[t-1,j]
            if(W[t,i,j]>vmax)
                vmax=W[t,i,j]
                g[t,i]=j
                V[t,i]=vmax
            end
        end
    end
end

# Plot results

#using Pkg
#Pkg.add("Plots")
using Plots, LaTeXStrings

# figure (1) : plot value function for some iterations
plot(K,V[1,:], label=L"V_{1}")
plot!(K,V[2,:], label=L"V_{2}")
plot!(K,V[3,:], label=L"V_{3}")
plot!(K,V[4,:], label=L"V_{4}")
plot!(K,V[8,:], label=L"V_{8}")
plot!(K,V[12,:], label=L"V_{12}")
plot!(K,V[160,:], label=L"V_{160}")

xlabel!("k")
ylabel!("V(k)")
title!("Value Function Iteration")


# figure (2) : plot optimal capital policy function
optimal_k_index = g[iter_max,:] # get index of optimal k_values
optimal_k_values = K[round.(Int, optimal_k_index)] # get corresponding value

plot(K, optimal_k_values, label=L"g(k)", legend=:bottomright)
plot!(K,K, label=L"45^o Line", color=:black)

xlabel!("k")
ylabel!("g(k)")
title!("Policy Function")

#= syntax notes
- the ! after plot! causes plots to be combined (multiple lines on one graph)
- the L"text" syntax creates a LaTeX field (used here L"V_{160}") =#