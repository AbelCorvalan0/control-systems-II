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
%% Simulate step response of transfer function.
%step(G)
%%% Conclusion.
% We've a inestable system due to two poles on origin.

%% Find the discrete open-loop transfer function Gd(s) of system with a
%% Zero-Order Retentor (Z0H) at the input. Implement sample time (Tm= 0.22 seg).
%Gd= minreal(c2d(G, Ts, 'zoh'))
Gd= c2d(G, Ts, 'zoh')
%% Poles and zeros map of continuous system.
%figure(2);
%pzmap(G);

%% Poles and zeros map of discrete system.
%disp('Zeros of the system: ')
%zero(Gd)
%disp('Poles of the system:')
%pole(Gd)
%% Zero-pole discrete system map.
%figure(3);
%pzmap(Gd)

%% Simulation with sample time Ts2= Ts*10.
%disp(sprintf('Discretized transfer function with Ts= %d', Ts2))
Gd1= c2d(G, Ts2, 'zoh')
%disp('Poles of the system Gd1: ')
%pole(Gd1)
%disp('Zeros of the system Gd1: ')
%zero(Gd1)
%% Poles and zeros map of the discrete system with
%% sample time Ts2= Ts*10.
%figure(4);
%pzmap(Gd1)

%% Obtain step response of discrete system. Determinate if this
%% is estable.
% Plot step response of continuous system.
figure(5);
step(Gd, 200)
% Plot step response of discrete system (Ts).
figure(6);
step(Gd1, 200)
% Plot step response of discrete system (10*Ts).
%figure(7);
%step(Gd1)

%%%%% Conclusion.
% The poles move towards the unit circle.

%%%%% CONCLUSIÓN
%% Sistema inestable ya que los polos se encuentran sobre
%% el círculo unitario en el plano Z.

%% Afecta el cero sobre el círculo unitario o en otras posiciones
%% a esta condición?

%% Por qué el hecho de los polos se encuentren dentro del círculo unitario
%% implica inestabilidad?

% Obtener la respuesta al escalon del sistema discreto y determinar
% si es estable.

%step(G)
%step(Gd, 11)
%step(Gd1, 11)

% Sistema discreto
% Determinar el tipo de sistema.

% Determinar la constante de error de posición Kp y el error ante un escalón y
% verificar mediante la respuesta al escalón de lazo cerrado del sistema
% discreto como se muestra.

%Kp= dcgain(Gd);
%Fd= feedback(Gd, 1)
%step(Fd, 20)

% Almaceno los valores de la respuesta al escalón
%[y, t]= step(Fd);
%s_response= [y, t];

% Obtengo los valores de amplitud de la respuesta al escalón.
%val_amplitud= s_response(:, 1);

% Realizo la diferencia de la ganancia
%e= zeros(length(val_amplitud), 1);

%for i= 1:1:length(val_amplitud)
%   e(i)= Kp-val_amplitud(i);
%end

%figure(4);
%plot(s_response(:, 2), e);
%title('Error entre ganancia DC (Kp) y respuesta al escalón');

% Verificar el error antes una rampa de entrada. ¿Converge o diverge.
% Explique la causa.

%t= 0:Tm:100*Tm;
%figure(5);
%lsim(Fd, t, t);
%ylim([0 2.5]);

%%%%% CONCLUSIÓN
%% Converge, Por qué? no sé.

% A lazo cerrado con realimentación unitaria

% Graficar el lugar de raíces del sistema continuo G(s) y del
% discreto Gd(s), indicando las ganancias criticas de estabilidad
% (si las hubiera)

%FTLC
%Gd= c2d(G, Tm, 'zoh');
%Gd1= c2d(G, 10*Tm, 'zoh');

% Lugar de raíces (Se obtiene con realimentación unitaria).
%rlocus(G);
%title(sprintf("Lugar de raíces G"));
%rlocus(Gd);
%title(sprintf("Lugar de raíces Gd con Tm= %.2f", Tm));

%disp('Estabilidad relativa si se aumenta 10 veces el tiempo de muestreo')
%title(sprintf("Lugar de raíces Gd con Tm= %.2f", 10*Tm));
%rlocus(Gd1);
