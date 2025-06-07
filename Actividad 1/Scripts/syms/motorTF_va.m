%% Sistemas de Control II - FCEFyN (UNC)

clear all; close all;
clc;

%% Define symbolic variables

syms Ra Laa Km Ki Bm J Tl
syms Ia wr va tita
syms s 

%% Model DC motor equations:

% Electrical characteristics.
% Armature coil can be described by:
% Laa: Inductance in series.
% Ra : Resistance in series.
% Va : Voltage source across the coil of the armature.

% Mechanical Characteristics.
% Ki: Torque constant.
% J : Intertia of the rotor and the equivalent mechanical load.
% Bm: Damping coefficient associated

eq1= s*Ia   == -(Ra/Laa)*Ia-(Km/Laa)*wr+(1/Laa)*va;
eq2= s*wr   == (Ki/J)*Ia-(Bm/J)*wr;  %-(1/J)*Tl;
eq3= s*tita == wr;

sol= solve([eq1, eq2, eq3],[wr, va, tita]);

disp('Original transfer function: ')
TFnn= sol.wr/sol.va
[num, den]= numden(TFnn)

%% Normalize Transfer Function
poly_num= Ki/(J*Laa);
poly_den= (((Bm*Ra+Ki*Km))/J*Laa)+((Bm*Laa+J*Ra)/J)*s+s^2;

N= coeffs(poly_num, s, 'All');
D= coeffs(poly_den, s, 'All');

disp('Normalized Transfer Function')
TF= poly_num/poly_den

%latex_str= latex(TF);
Ki= 0.1448
Bm= 0
Ra= 2.4134

eq4= J*Laa== 0.0006374
eq5= (Bm*Laa+J*Ra)== 0.09658
eq6= (Bm*Ra+Ki*Km)== 1