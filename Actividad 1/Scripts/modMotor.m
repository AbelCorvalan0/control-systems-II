function [X] = modMotor(t_etapa, xant, accion, Tl)
% modMotor (función obtenida de los apuntes de clase)
%   Esta función modela de forma lineal el funcionamiento de un motor de CC
%   a partir de las ecuaciones
%   1) ia_p   = -Ra/Laa*ia - Km/Laa*wr + 1/Laa*Va
%   2) wr_p   = Ki/J*ia - Bm/J*wr - 1/J*TL
%   3) tita_p = wr
%
%   Args:
%   - t_etapa: duración de la etapa actual
%   - xant   : valores anteriores de las salidas de interés
%   - accion : entrada de control sobre el motor
%
%   Output:
%   - X: valores de las salidas de interés en el tiempo actual

% Parámetros del motor a modelar
% Laa  = 0.002511;
% Ki   = 3.765;
% J    = 0.02056;
% Km   = 0.2485;
% B    = 0.026;
% Ra   = 2.415;
Laa = 0.0047;       % [Hy].
J   = 0.00233;      % [kg.m^2].
Ra  = 2.27;         % [Ohm].
B   = 0.00131;      % [N.m/(rad/seg)].
Ki  = 0.25;         % 
Km  = 0.25;         % 

Va = accion;          % tensión de armadura
h  = 1e-7;            % paso de integración

TL = Tl;

% Asignación de los valores actuales de las variables de interés para
% iniciar el proceso de cálculo
wr    = xant(1);
wrp   = xant(2); 
ia    = xant(3);
theta = xant(4);

for ii=1:t_etapa/h
    % Cálculo de las derivadas
    wrpp = (-wrp*(Ra*J + Laa*B) - wr*(Ra*B + Ki*Km) + Va*Ki)/(J*Laa);
    iap  = (-Ra*ia - Km*wr + Va)/Laa;

    % Aplicación del método de Euler para la estimación de los valores de
    % las variables de interés en el tiempo posterior
    wrp  = wrp + h*wrpp;
    wrp  = wrp - ((1/J)*TL);          % reducción por acción del torque
    ia   = ia + iap*h;                % cálculo de la corriente de armadura
    wr   = wr + h*wrp;                % cálculo de la velocidad angular
    theta = theta + h*wr;             % cálculo de la posición 
end

% Vector de salida
X = [wr, wrp, ia, theta];
end