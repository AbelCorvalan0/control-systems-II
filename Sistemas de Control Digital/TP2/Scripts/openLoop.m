% clear all; close all;
% clc

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

% Kfin= Kp + Kd
% c= Kd/(Kp + Kd)= Kd /K

%% Cálculo de ganacias controlador PD.

K= 0.23357
c= 0.7663
Kd= c*K
Kp= K - Kd

% Paleta de colores personalizada
col_salida      = [0 0.4470 0.7410];  % azul
col_controlador = [0.4660 0.6740 0.1880]; % verde
col_error       = [0.8500 0.3250 0.0980];  % rojo
col_derivativa  = [0.4940 0.1840 0.5560];  % violeta
col_integral    = [0.9290 0.6940 0.1250];  % mostaza
col_proporcional= [0.3010 0.7450 0.9330];  % celeste

subplot(3,2,1);
plot(tout, yout(:,1), 'Color', col_salida, 'LineWidth', 1.5);
title('1. Salida del sistema');
xlabel('Tiempo [s]'); ylabel('Salida');
grid on;

subplot(3,2,2);
plot(tout, yout(:,2), 'Color', col_controlador, 'LineWidth', 1.5);
title('2. Salida del controlador');
xlabel('Tiempo [s]'); ylabel('u(t)');
grid on;

subplot(3,2,3);
plot(tout, yout(:,3), 'Color', col_error, 'LineWidth', 1.5);
title('3. Error');
xlabel('Tiempo [s]'); ylabel('e(t)');
grid on;

subplot(3,2,4);
plot(tout, yout(:,4), 'Color', col_derivativa, 'LineWidth', 1.5);
title('4. Acción derivativa');
xlabel('Tiempo [s]'); ylabel('D(t)');
grid on;
 
subplot(3,2,5);
plot(tout, yout(:,5), 'Color', col_integral, 'LineWidth', 1.5);
title('5. Acción integral');
xlabel('Tiempo [s]'); ylabel('I(t)');
grid on;

subplot(3,2,6);
plot(tout, yout(:,6), 'Color', col_proporcional, 'LineWidth', 1.5);
title('6. Acción proporcional');
xlabel('Tiempo [s]'); ylabel('P(t)');
grid on;