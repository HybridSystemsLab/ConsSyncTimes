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
H = -.4;                                           % Controller Flow Gain
v = [.5,1.5];                                      % Communication Interval
G = [0 1 1 0 1; 1 0 1 0 0; 1 0 0 1 0; ....
    0 0 1 0 1; 1 0 1 1 0];                         % Adjacency matrix
L = diag(sum(G')) - G;                             % Laplacian

% Initial Conditions (ICs)
x0 = [1 -1 2 -2 0]';                               % State ICs
eta0 = [0 -3 1 -4 -1]';                            % Hybrid Controller ICs
Tau0 = .2;                                         % Event Timer IC
X0 = [x0; eta0;Tau0];

% Simulation horizon
TSPAN=[0 20];                                      % Flow Horizon
JSPAN = [0 1000];                                  % Jump Horizon

%Simulation Options
options = odeset('RelTol',1e-6,'MaxStep',.1);

% Simulate
[t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,1,options,1);

% Found using a convex optimization package with the above paramters (cvx)
P1 = [25.2540246909757,-8.46826567736576e-16,9.50183823241938e-16,-6.88601530807639e-16,3.16134682382400,-3.24640445205895e-16,-7.55065294773646e-16,-2.44110281446369e-16;-8.46826567736576e-16,19.3480884572871,-2.93703663634259e-15,-4.54323491408958e-16,2.66521886762185e-16,3.86075517294809,6.11578082326176e-16,1.79364996230616e-16;9.50183823241938e-16,-2.93703663634259e-15,17.4537714780361,-3.10291543798892e-15,3.87443550730015e-16,-1.67364000843291e-15,2.61992199046472,-1.78094868531848e-15;-6.88601530807639e-16,-4.54323491408958e-16,-3.10291543798892e-15,19.3480884575222,8.46722104827283e-17,-4.43406287992792e-16,5.48142011262904e-16,3.86075517302153;3.16134682382400,2.66521886762185e-16,3.87443550730015e-16,8.46722104827283e-17,5.84371626033440,3.75585786690253e-16,8.35399701878592e-16,8.69060209232409e-17;-3.24640445205895e-16,3.86075517294809,-1.67364000843291e-15,-4.43406287992792e-16,3.75585786690253e-16,8.74869301825895,2.59713301987607e-16,-2.33636276047575e-16;-7.55065294773646e-16,6.11578082326176e-16,2.61992199046472,5.48142011262904e-16,8.35399701878592e-16,2.59713301987607e-16,11.2704892302197,4.26433352315976e-16;-2.44110281446369e-16,1.79364996230616e-16,-1.78094868531848e-15,3.86075517302153,8.69060209232409e-17,-2.33636276047575e-16,4.26433352315976e-16,8.74869301850018];


% Evaluating the Lyapunov function
bary = NaN(length(y),2*N);
V = NaN(length(y),1);
expmAf = @(tau) expm(Af*tau);
expmAft = @(tau) expm(Af'*tau);
BTi = kron(eye(2),T^-1);
for i = 1:length(t)
    bary(i,:) = BTi * y(i,1:end-1)';
    V(i) = bary(i,[2:N, N+2:end])*expmAft(y(i,end))*P1*expmAf(y(i,end))*bary(i,[2:N, N+2:end])';
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

