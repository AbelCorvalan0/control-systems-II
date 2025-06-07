clear all; close all;
clc

%%%%% Import libraries.
pkg load control
pkg load symbolic
syms s real

%%%%% Define values.
%% Define pole values.
p1 = 0;
p2 = 0;
%% Define zero value.
z1 = -10;
%% Define gain value.
K  = 5;
% (override) = 5
% 2% Time= 3
% error  = 0;
%% Define sample time.
Ts = 0.22;
Ts2= Ts*10;
%%%%% System in open-loop condition.
%% Continuous Transfer Function.
disp('Transfer Function: ')
G= minreal(zpk([z1], [p1 p2], K))

Gd= c2d(G, Ts, 'zoh')

kp= dcgain(Gd)
F= feedback(Gd, 1)
step(F); close;

%% Verifying of error due to ramp input.
% Build ramp input.
%t= 0:Ts:100*Ts;
% Ramp response simulation of close-loop system.
%ramp(F)