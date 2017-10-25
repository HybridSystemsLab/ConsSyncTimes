%{  
    File:    run4_consensus.m
    Author:  Sean Phillips
    
    Utility: Main run file
    
    Project: Hybrid simulation of static consensus of 
             multiple agent network systems with 
             global communication events.
    
    Dependencies:
    
    C.m
    D.m
    f.m
    g.m
    hybridsolver.m
    plotflows.m

%}

global K v N G L

% System Parameters
K = .3;                                            % Controller Jump Gain
H = 0;                                           % Controller Flow Gain
v = [.5,1.5];                                      % Communication Interval
G = [0 1 1 0 1; 1 0 1 0 0; 1 0 0 1 0; ....
    0 0 1 0 1; 1 0 1 1 0];                         % Adjacency matrix
L = diag(sum(G')) - G;                             % Laplacian
N = size(G,2);                                     % Number of Agents

% Initial Conditions (ICs)
x0 = [1 -1 2 -2 0]';                               % State ICs
eta0 = [0 -3 1 -4 -1]';                            % Hybrid Controller ICs
Tau0 = .2;                                         % Event Timer IC
X0 = [x0; eta0;Tau0];                              % Initial Conditions

% Simulation horizon
TSPAN=[0 20];                                      % Flow Horizon
JSPAN = [0 1000];                                  % Jump Horizon

%Simulation Options
options = odeset('RelTol',1e-6,'MaxStep',.1);

% Simulate
[t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,1,options,1);

% Found using a convex optimization package with the above paramters (cvx)
P1 =   [    33.6103    0.0000   -0.0000    0.0000    4.2030   -0.0000    0.0000   -0.0000;...
    0.0000   28.6082    0.0000    0.0000   -0.0000    5.7271   -0.0000    0.0000;...
   -0.0000    0.0000   25.3450   -0.0000    0.0000    0.0000    4.7508   -0.0000;...
    0.0000    0.0000   -0.0000   28.6082    0.0000   -0.0000    0.0000    5.7271;...
    4.2030   -0.0000    0.0000    0.0000    7.0216   -0.0000    0.0000    0.0000;...
   -0.0000    5.7271    0.0000   -0.0000   -0.0000   11.1316    0.0000   -0.0000;...
    0.0000   -0.0000    4.7508    0.0000    0.0000    0.0000   14.9554   -0.0000;...
   -0.0000    0.0000   -0.0000    5.7271    0.0000   -0.0000   -0.0000   11.1316];

% Jordan Matrix form of Laplacian
T = [    0.4472   -0.6914    0.3015    0.4157    0.6040;...
    0.4472    0.1058   -0.6030    0.7763   -0.1060;...
    0.4472    0.4360    0.3015   -0.0942   -0.4980;...
    0.4472   -0.3612   -0.6030   -0.4547   -0.1060;...
    0.4472    0.4360    0.3015   -0.0942    0.6040];

% Evaluating the Lyapunov function
bary = NaN(length(y),2*N);                              % Initialize bar y
V = NaN(length(y),1);                                   % Initialize V
Z = zeros(N-1);                                         % Zero matrix size N
I = eye(N-1);                                           % Identity matrix size N
Af = [Z, I; Z, Z];                                      % Flow matrix

BTi = kron(eye(2),T);
for i = 1:length(t)                                     % Evaluate V over solution
    bary(i,:) = BTi * y(i,1:end-1)';
    baryi = [bary(i,2:N), bary(i,N+2:2*N)]';
    tau = y(i,end);
    V(i) = baryi'*expm(Af*tau)'*P1*expm(Af*tau)*baryi;
end

% Plotting the states, the timer and the Lyapunov function over the
% solutions
figure(1)
clf
subplot(311)
h = plot(t,y(:,1:N),'LineWidth',2);
hold on
legend('', '', '', '', '')
axis([0,5,-3,3])
subplot(312)
plotflows(t,j,y(:,11),[0 0.4470 0.7410])
subplot(313)
plotflows(t,j,V,[0.4940 0.1840 0.5560])

