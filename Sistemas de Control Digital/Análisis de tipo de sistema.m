%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Tarea 1 - Análisis de tipo de sistema %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

pkg load symbolic
syms z real

%%% Función transferencia en el tiempo discreto.
Gd= (0.01852*z+0.01693)/((z^2)-1.749*z+0.7634)
%Gd= (0.3933*z-0.3933)/((z^2)+1.749*z+0.7634)

% Tiempo de muestreo
Tm= 0.09

% Defino el integrador
%it= (1/(z-1))
it= (z-1)

% Determino el tipo de sistema.
%Kp= lim Gd(z)
%    z->1
Kp= Gd

%Kv= lim (1/(z-1))*G(s)
%    z->1
Kv= simplify((it)*Gd)

%Ka= lim ((1/((z-1)^2)))*G(s)
%    z->1
Ka= simplify((it^2)*Gd)

%evaluo z= 1.
z= 1
disp("El valor de Kp es:")
eval(Kp)
disp("El valor de Kv es:")
eval(Kv)
disp("El valor de Ka es:")
eval(Ka)

%El error en estado estable para una entrada en escalon es ess=0=1/(1+Kp)
%El error en estado estable para una entrada rampa es ess=cte=1/Kv
%El error en estado estable para una entrada parabólica es ess=inf=1/Ka

ess_e= eval(1/(1+Kp))
ess_r= eval((1/Kv))
ess_p= eval((1/Ka))

disp("El sistema es de tipo 0")
