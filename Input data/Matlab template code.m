%  This matlab program performs the value function iteration for 
%  social planner's problem in the neo-classical growth model.
%  Serdar Ozkan, ECO 2061, Winter 2016, serdar.ozkan@utoronto.ca
%%
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
