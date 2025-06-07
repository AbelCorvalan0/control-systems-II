function [X] = timesCalculator(polM, polm)

%% Integration time calculating.
tR   = log(0.95)/(polM);
tint = tR/3;

%% Simulation time calculating.
tL= log(0.05)/polm;
%% Output vector.
X = [tint, tsim];
end