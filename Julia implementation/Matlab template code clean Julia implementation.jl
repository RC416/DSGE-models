
iter_max=160; # Choose your maximum number of iteration to make sure that V converges

α=0.3;
β=0.65;
k_bar=(α*β)^(1/(1-α))
# Make sure that your grid points for k include the steady state value of k
K=0.05:0.025:0.15; 
# K=0.05:0.01:0.50; # This is a finer grid points.
β
N = length(K)

# initialize value function 
V = zeros(iter_max, N)
V[1,:]=zeros(1,N); # This is my initial guess. 

# initialize capital policy function and temporary value function for the loop (W)
g = zeros(iter_max, N)
W = zeros(iter_max, N, N)

for t=2:iter_max
    for i=1:N
        vmax=-100000000;
        for j=1:N
            W[t,i,j]=log(K[i]^α-K[j])+β*V[t-1,j];
            if(W[t,i,j]>vmax)
                vmax=W[t,i,j];
                g[t,i]=j;
                V[t,i]=vmax;
            end
        end
    end
end

