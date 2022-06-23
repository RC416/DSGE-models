% Base Model implementation in Matlab

% Assign parameter values.
alpha = 0.400;
beta  = 0.987;
delta = 1.000;
number_of_iterations = 1000;

% Calculate the steady-state level of capital.
k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ^ (1/(alpha-1));

% Create a range of capital values around steady-state (+/- 50%).
number_of_k_values = 201;
k_low_pct = 0.50;
k_high_pct = 1.50;
k_values = linspace(k_low_pct*k_steady, k_high_pct*k_steady, number_of_k_values);

% Initialize Value Function and Policy Function (as arrays).
Value_Function = zeros(number_of_iterations, number_of_k_values);
Policy_Function = zeros(number_of_iterations, number_of_k_values);


% Perform value function iteration.
for iteration = 2:(number_of_iterations)
    
    for ind_kt0 = 1:number_of_k_values      % for each level of starting captial...

        v_max = -inf; % placeholder for candidate value function value

        for ind_kt1 = 1:number_of_k_values  % ...check all next period capital choices

            % Calculate value of Value Function for given starting capital and next period capital choice.
            New_Value_Function_Value = log((k_values(ind_kt0)^alpha) - k_values(ind_kt1) + (1-delta)*k_values(ind_kt0))...
                + beta*Value_Function(iteration-1, ind_kt1);
            
            % Check if this capital choice gives highest Value Function value.
            if New_Value_Function_Value > v_max
                
                % Update value of v_max, Value Function, and Policy Function.
                v_max = New_Value_Function_Value;
                Value_Function(iteration, ind_kt0) = v_max;     % store value
                Policy_Function(iteration, ind_kt0) = ind_kt1;  % store index
            end
        end
    end
end


% Plot various iterations of the Value Function.


% Plot final Policy Function.















clear
clc
iter_max=160; % Choose your maximum number of iteration to make sure that V converges
alpha=0.3;
beta=0.65;
k_bar=(alpha*beta)^(1/(1-alpha))
% Make sure that your grid points for k include the steady state value of k
K=[0.05:0.025:0.15]; 
% K=[0.05:0.01:0.50]; % This is a finer grid points.
[m,N]=size(K);
V(1,:)=zeros(1,N); % This is my initial guess. 

for t=2:iter_max
    for i=1:N
        vmax=-100000000;
        for j=1:N
            W(t,i,j)=log(K(i)^alpha-K(j))+beta*V(t-1,j);
            if(W(t,i,j)>vmax)
                vmax=W(t,i,j);
                g(t,i)=j; % Policy function
                V(t,i)=vmax; % Value function
            end
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following are a bunch of defaults for plots that make your figures
% look nice.
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultLineMarkerSize',8);
set(0,'defaultlinelinewidth',4);
set(0,'defaultTextFontSize',20);
set(0,'DefaultAxesFontSize',24);
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultAxesFontName','Times New Roman');
lfsize=24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1) % Plotting value functions for some iterations
hold off
plot(K,V(1,:),'+-')
hold on
plot(K,V(2,:),'+-')
plot(K,V(3,:),'+-')
plot(K,V(4,:),'+-')
plot(K,V(8,:),'+-')
plot(K,V(12,:),'+-')
plot(K,V(160,:),'+-')
grid on
xlabel('$k$')
ylabel('$V(k)$')
title('Value Function Iteration')
leg=legend('V_0','V_1','V_2','V_4','V_8','V_{12}','V_N');
set(leg,'FontSize',lfsize);
axis tight
%%
figure(2) % Plotting policy function. 
hold off
plot(K,K(g(iter_max,:)),'+-')
hold on
plot(K,K)
grid on
xlabel('$k$')
ylabel('$g(k)$')
title('Policy Function')
leg=legend('g(k)','45^o Line');
set(leg,'FontSize',lfsize);
axis tight
