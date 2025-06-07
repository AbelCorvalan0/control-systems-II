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

%%
Gd1= c2d(G, Ts2, 'zoh')

%%
figure(1);
rlocus(G)

% figure(2);
% rlocus(Gd)

% figure(3);
% rlocus(Gd1)

%% Calculate gain margin of G(s).
[Gm, ~, ~]= margin(Gd1);
Kcrit= Gm;    % Critical K for unit feedback.
disp(["Critical gain K: ", num2str(Kcrit)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Relative stability analysis.

% 2. Análisis de estabilidad relativa (márgenes de ganancia y fase)
[GM_linear, PM, Wcg, Wcp] = margin(Gd1); % GM en escala lineal
GM_dB = 20 * log10(GM_linear);           % Convertir GM a dB

disp(['Margen de Ganancia (GM) en dB: ', num2str(GM_dB), ' dB']);
disp(['Margen de Fase (PM): ', num2str(PM), ' grados']);
disp(['Frecuencia de cruce de ganancia (Wcg): ', num2str(Wcg), ' rad/s']);
disp(['Frecuencia de cruce de fase (Wcp): ', num2str(Wcp), ' rad/s']);

% Interpretación de GM en dB:
% GM > 0 dB: El sistema puede tolerar un aumento de ganancia antes de volverse inestable.
% GM = 0 dB: Límite de estabilidad.
% GM < 0 dB: El sistema es inestable para ganancias nominales.

disp('El sistema es inestable para ganancias nominales.')

% 3. Diagramas de Bode para visualización (GM y PM)
figure;
bode(Gd1);
%grid on;
title('Diagrama de Bode del Sistema Discreto');

% 4. Mapa de polos y ceros
figure;
pzmap(Gd);
figure;
pzmap(Gd1);
%grid on;
title('Mapa de Polos y Ceros del Sistema Discreto');