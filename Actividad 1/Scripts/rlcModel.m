function [X]= rlcModel(r, l, c, u)
clear all; close all; clc;

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

%% Operation point.
% Initial values x(0)
Il(1)  = 0;        % Current.
Vrl(1) = 0;        % Voltage through R. ??
y(1)   = 0;        % Output.

Xop    = [0, 0]';
x      = [Il(1) Vrl(1)]'; 

R= r; L= l; C= c;

%% State Space.
A= [-R/L -1/L; 1/C 0]   % State matrix.
B= [1/L; 0]             % Input matrix.
c= [R 0]                % Output matrix.
D= 0;                    % Direct transmission matrix.

%% Obtain Transfer Function.
% [num, den]= ss2tf(A, B, c, D);
% G= tf(num, den)     % Build system transfer function. 
% p= pole(G)
% % Fastest dynamic pole.
% p(1)
% % Slower dynamic pole.
% p(2)

% tint = log(0.95)/p(1);
% tsim = log(0.05)/p(2);
% T= 20e-3; 

xp= A*(x-Xop)+B*u;
x = x+ xp*tint;
Y= c*(x-Xop)+ D*u;

ii= i+1;
y(ii)   = Y(1);
Il(ii)  = x(1);
Vrl(ii) = x(2);

X= [y];
end
