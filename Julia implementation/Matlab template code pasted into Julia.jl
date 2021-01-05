#=

Template code from Matlab pasted into Julia

It runs with only the following changes:
 1. changed comment-out character from % to #
 2. removed brackets around initiation of the list K (line 7/8)
 3. changed initiation of m,N (line 11/12)
 4. changed subscripting bracks to square: V(1,1) to V[1,1] (lines 16, 25, 26-29)
 5. pre-initiated all policy funciton V, g, W (lines 15, 18, 19)
 
=#

iter_max=160; # Choose your maximum number of iteration to make sure that V converges

alpha=0.3;
beta=0.65;
k_bar=(alpha*beta)^(1/(1-alpha))
# Make sure that your grid points for k include the steady state value of k
# K=[0.05:0.025:0.15]; 
K=0.05:0.025:0.15; 
# K=[0.05:0.01:0.50]; # This is a finer grid points.

# [m,N]=size(K);
N = length(K)

# V(1,:)=zeros(1,N); # This is my initial guess. 
V = zeros(iter_max, N)
V[1,:]=zeros(1,N); # This is my initial guess. 

g = zeros(iter_max, N)
W = zeros(iter_max, N, N)

for t=2:iter_max
    for i=1:N
        vmax=-100000000;
        for j=1:N
            W[t,i,j]=log(K[i]^alpha-K[j])+beta*V[t-1,j];    # changed curly brackets to square for subscripting
            if(W[t,i,j]>vmax)
                vmax=W[t,i,j];
                g[t,i]=j; # Policy function
                V[t,i]=vmax; # Value function
            end
        end
    end
end


using Pkg
Pkg.add("plots")
using Plots

plot(K,V[1,:],"+-")