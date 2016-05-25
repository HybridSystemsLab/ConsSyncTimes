function xplus = g(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               Author: Torstein Ingebrigtsen Bø
%
% Project: Simulation of a hybrid system (bouncing ball)
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K v N Adj m L noisy

xg = x(1:N);
% Assume All-to-All Connection

    tauplus =v(1) + (v(2)-v(1))*rand(1);
if noisy == 1
    m = .1*rand(N,1);
    etaplus =-K*L*xg - K*L*m;
else
    etaplus =-K*L*xg;
end

xplus = [x(1:N);etaplus;tauplus];

end