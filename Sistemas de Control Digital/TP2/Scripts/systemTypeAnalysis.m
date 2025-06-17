%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Tarea 1 - Análisis de tipo de sistema %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

pkg load symbolic
syms z real

%%% Función transferencia en el tiempo discreto.
Gd= (2.31*z+0.11)/((z^2)-2*z+1)
%Gd= (0.3933*z-0.3933)/((z^2)+1.749*z+0.7634)

% Tiempo de muestreo
Ts= 0.22

% Defino el integrador
%it= (1/(z-1))
it= (z-1)

% Determino el tipo de sistema.
%Kp= lim Gd(z)
%    z->1
Kp= Gd

%Kv= lim ((z-1)/Ts)*G(s)
%    z->1
Kv= simplify(((it)/(Ts))*Gd)

%Ka= lim (((z-1)^2/(Ts^(2))))*G(s)
%    z->1
Ka= simplify((((it)^(2))/(Ts^2))*Gd)

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

disp("El sistema es de tipo 2")

%Kp= cte sistema es de tipo 0
%Kp= inf sistema de tipo 1 o mayor.

%Kv= cte sistema de tipo 1.
%Kv= 0   sistema de tipo 0.
%Kv= inf sistema de tipo 2 o mayor.

%Ka= cte sistema de tipo 2.
%Ka= 0   sistema tipo menor a 2.