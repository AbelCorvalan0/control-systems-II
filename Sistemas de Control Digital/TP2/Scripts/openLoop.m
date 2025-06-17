clear all; close all;
clc

%%%%% Import libraries.
%pkg load control
%pkg load symbolic
%syms s real

% Ganancia  = 5;
% Sobrepaso = 5;
% tiempo 2% = 3;
% ts        = 0.22;

%%%%% Define values.
%% Define pole values.
p1 = 0;
p2 = 0;
%% Define zero value.
z1 = -10;
%% Define gain value.
K  = 5;
% (override) = 5.
% 2% Time= 3.
% error  = 0;
%% Define sample time.
Ts = 0.22;
Ts2= Ts*10;
%%%%% System in open-loop condition.
%% Continuous Transfer Function.
disp('Transfer Function: ')
G= minreal(zpk([z1], [p1 p2], K))
%% Simulate step response of transfer function.
%step(G)
%%% Conclusion.
% We've a inestable system due to two poles on origin.

%% Discretized system.
Gd= c2d(G, Ts, 'zoh')

step(Gd)
pole(Gd)
zero(Gd)

%% Calculation of psita, natural frequency (\omega_{0})
% Overshoot.
S= 5;
psita= -log(S/100)/(sqrt(pi^(2)+(log(S/100))^(2)))
% psita= 0.6901 
% 0<psita<1 decrement oscillations.|

% tR(2%).
tR = 3;
w0 = 4/(psita*tR)

wd= w0*sqrt(1-(psita)^(2))

td= (2*pi)/(wd)

%% Cálculo de muestras por ciclo de la frecuencia amortiguada wd.
Ts= 0.22;
m= td/Ts

%% Mediante la equivalencia de planos s y z
% determinar la ubicación de los polos deseados en el plano z.

r= exp(-psita*w0*Ts)

ang= wd*Ts

%% Build close-loop system.
% SISOTOOL.
%sisotool(Gd)
%kp= dcgain(Gd)

F= feedback(Gd, 1)

%C
%F= feedback(C*Gd, 1)
pole(F)
zero(F)
pzmap(F)
step(F); 
%close;

%% Conclusions

% Lead compensator: doesn't verify maximium overshoot 
% (takes 1.35 amplitude and allowed value is 1.05).
% Test cancelling pole between marks.

