clear all; close all; clc;

Laa= 366e-6;    % [Hy].
J  = 5e-9;      % [kg.m^2].
Ra = 55.6;      % [Ohm].
Bm = 0;         % [N.m/(rad/seg)].
Ki = 6.49e-3;   % 
Km = 6.53e-3;   % 

s   = tf('s');
num = [Ki/(J*Laa)];
den = [1,((Bm*Laa + J*Ra))/J, (Laa*(Bm*Ra + Ki*Km))/J];

TF= tf(num, den)

[A, B, C, D]= tf2ss(num, den);
% Obtain poles of Transfer Function.
p= pole(TF);
%rlocus(TF);

%% Simulation

h= 10e-7;

% Simulation time calculate.
tsim= log(0.05)/p(2);
ts= round(tsim)

c2p= 4;   % cycles2plot
t= 0: h: c2p*ts;
u= linspace(0, length(t), h);
size(t)
cycles= 0;
point = 0;
for i=1: (c2p*ts)
    if point>ts
        u(i)= 12;
    end
    point= point+1;
    % Xp= A*(x-xop)+ B*(u);
    % y= C*x+D*u;
end
size(u)
% figure(1)
% plot(t, u)
% grid on;