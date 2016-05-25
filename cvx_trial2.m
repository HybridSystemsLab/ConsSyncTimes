
% close all
clear variables


tic;

global K v N Adj A B n L noisy
noisy = 0;
A = 0;
B = 1;
n = length(A);

N = 5;

Adj = [0 1 1 0 1 ;
        1 0 1 0 0 ;
        1 0 0 1 0 ;
        0 0 1 0 1 ; 
        1 0 1 1 0];

% Adj = [0 1 1 0 1 ;
%         1 0 1 0 0 ;
%         1 0 0 1 0 ;
%         0 0 1 0 1 ; 
%         1 0 1 1 0];
% Adj = [0 1 1 1 1; 
%        1 0 1 1 1; 
%        1 1 0 1 1; 
%        1 1 1 0 1; 
%        1 1 1 1 0]; 
% Adj = [0 1 0 0 0; 
%        1 0 1 0 1; 
%        0 1 0 1 0; 
%        0 0 1 0 0; 
%        0 1 0 0 0];
K = .3;

v = [0.1,1];
% beta = 0.15;

Din = diag(sum(Adj'));
L = Din - Adj;

[T,Q] = eig(L,'vector');
indx = find(abs(Q) < .0001);
T(:,[1,indx]) = T(:,[indx,1]); % If indx neq 0, then it will switch columns

F = T^-1*L*T;

barL = F(2:end,2:end);

Z = zeros(length(barL));
I = eye(length(barL));

Af = [Z I; Z Z];
Ag = [I Z; -K*barL Z];


A1low = expm(Af*v(1))*Ag;
A1high = expm(Af*v(2))*Ag;
eigA1high = max(abs(eig(A1high)))
eigA1low = max(abs(eig(A1low)))

if eigA1high > 1
    error('eigA1high > 1')
end

cvx_begin sdp quiet
    variable P1(size(Ag)) symmetric 
    minimize (0)
    subject to
    A1high'*P1*A1high - P1 <= 0;
    A1low'*P1*A1low - P1 <= 0;
    P1 >= eye(size(Ag));
cvx_end


% simulation horizon
TSPAN=[0 5];
JSPAN = [0 1000];

% rule for jumps
% rule = 1 -> priority for jumps`
% rule = 2 -> priority for flows
rule = 1;



options = odeset('RelTol',1e-2,'MaxStep',.1);

% if isnan(P1(1,1))
%     figure(1)
%     subplot(211)
%     plot(t, (eye(N) - Adj)*y(:,1:N)')
%     subplot(212)
%     plot(t,y(:,1:N))
%     error('Could not find P')
%     
% else
%     display('Good P')
%     display(abs(eig(P1)))
% end

%% Evaluate Lyapunov Function
% 


x0 = [1 -1 2 -2 0]';
eta0 = [0 -3 1 -4 -1]';
Tau0 = .2;

X0 = [x0; eta0;Tau0];

[t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,rule,options);

bary = NaN(length(y),2*N);
V = NaN(length(y),1);
expmAf = @(tau) expm(Af*tau);
expmAft = @(tau) expm(Af'*tau);
BTi = kron(eye(2),T^-1);
for i = 1:length(t)
    bary(i,:) = BTi * y(i,1:end-1)';
    V(i) = bary(i,[2:N, N+2:end])*expmAft(y(i,end))*P1*expmAf(y(i,end))*bary(i,[2:N, N+2:end])';
end

%% Evaluate local error

bary = NaN(length(y),2*N);
V = NaN(length(y),1);
BTi = kron(eye(2),T^-1);
for i = 1:length(t)
    bary(i,:) = BTi * y(i,1:end-1)';
    V(i) = bary(i,[2:N, N+2:end])*expmAft(y(i,end))*P1*expmAf(y(i,end))*bary(i,[2:N, N+2:end])';
end


%% Plot For Example 1
clf

subplot(311)
h = plot(t,y(:,1:N),'LineWidth',2)
hold on
legend('', '', '', '', '')
axis([0,5,-3,3])
subplot(312)
plotflows(t,j,y(:,11),[0 0.4470 0.7410])
subplot(313)
plotflows(t,j,V,[0.4940 0.1840 0.5560])
% axis([0,5,-5,5])

%% Plot for Example 2
% clf
% avg = 1/N*(ones(1,N)*y(1,1:5)' + ones(1,N)*y(1,6:10)'*y(1,11));
% 
% % subplot(211)
% plot(t,y(:,1:5),'LineWidth',2)
% hold on
% plot(t,avg*ones(size(t)),'--k','LineWidth',2)
% legend('', '', '', '', '','')
% subplot(212)
% plot(t,log10(V))
% axis([0,10,0,1100])

%% Plot for Example 3 - 2 agent case

for ii = 1:10
    [t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,rule,options);

    avg = 1/N*(ones(1,N)*y(1,1:2)' + ones(1,N)*y(1,3:4)'*y(1,5));

    figure(1)
    subplot(211)
    plot(t,y(:,1),'LineWidth',2)
    hold on
    plot(t,avg*ones(size(t)),'--k','LineWidth',2)
    subplot(212)
    plot(t,y(:,2),'LineWidth',2)
    hold on
    plot(t,avg*ones(size(t)),'--k','LineWidth',2)
    drawnow
end

%% Plot for partial Point-wise

% figure(1)
% plot([0,.8],[0,.8],'LineWidth',2)
% figure(2)
% figure(2)
% subplot(411)
% hold on 
% plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*.6,'k--','LineWidth',2)
% plot(t,ones(size(t))*.4,'k','LineWidth',2)
% subplot(412)
% hold on 
% plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*.6,'k--','LineWidth',2)
% plot(t,ones(size(t))*.4,'k','LineWidth',2)
% subplot(413)
% hold on 
% plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*-.2,'k--','LineWidth',2)
% plot(t,ones(size(t))*0,'k','LineWidth',2)
% subplot(414)
% hold on 
% plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*-.2,'k--','LineWidth',2)
% plot(t,ones(size(t))*0,'k','LineWidth',2)


for ii = 1:25
    xrand = .2 - .4*rand(2,1);
    Tau0 = rand(1);
    x0 = [.4, .4]' + xrand;
    eta0 = .2 - .4*rand(2,1);
    X0 = [x0;eta0;Tau0];
   
    [t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,rule,options);

    figure(2)
    subplot(411)
    hold on 
    plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*.6,'k--','LineWidth',2)
    plot(t,ones(size(t))*.4,'k','LineWidth',2)
    subplot(412)
    hold on 
    plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*.6,'k--','LineWidth',2)
    plot(t,ones(size(t))*.4,'k','LineWidth',2)
    subplot(413)
    hold on 
    plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*-.2,'k--','LineWidth',2)
    plot(t,ones(size(t))*0,'k','LineWidth',2)
    subplot(414)
    hold on 
    plot(t,ones(size(t))*.2,'k--',t,ones(size(t))*-.2,'k--','LineWidth',2)
    plot(t,ones(size(t))*0,'k','LineWidth',2)
    
    avg = 1/N*(ones(1,N)*y(1,1:2)' + ones(1,N)*y(1,3:4)'*y(1,5));
%     figure(1)
%     hold on
%     plot(y(:,1),y(:,2),'LineWidth',2)
    figure(2)
    subplot(411)
    hold on 
    plot(t,y(:,1),'LineWidth',2)
    subplot(412)
    hold on 
    plot(t,y(:,2),'LineWidth',2)
    subplot(413)
    hold on 
    plot(t,y(:,3),'LineWidth',2)
    subplot(414)
    hold on 
    plot(t,y(:,4),'LineWidth',2)
    drawnow;
end

expmAfPexpmAf = @(tau) expm(Af'*tau)*P1*expm(Af*tau);

eigmin = 100;
eigmax = 0;
for i = 0:.05:v(2)
    if min(eig(expmAfPexpmAf(i))) < eigmin
        eigmin = min(abs(eig(expmAfPexpmAf(i))));
    end
    if max(eig(expmAfPexpmAf(i))) > eigmax
        eigmax = max(abs(eig(expmAfPexpmAf(i))));
    end
end


e = sqrt(eigmax/eigmin)/4;
subplot(411)
plot(t,ones(size(t))*.4-e,'r--',t,ones(size(t))*.4+e,'r--','LineWidth',2)
subplot(412)
hold on 
plot(t,ones(size(t))*.4-e,'r--',t,ones(size(t))*.4+e,'r--','LineWidth',2)
subplot(413)
hold on 
plot(t,ones(size(t))*-e,'r--',t,ones(size(t))*e,'r--','LineWidth',2)
subplot(414)
hold on 
plot(t,ones(size(t))*-e,'r--',t,ones(size(t))*e,'r--','LineWidth',2)




