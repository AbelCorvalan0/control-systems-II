clear all; close all;
clc

pkg load control
pkg load symbolic
syms s real

%Función transferencia continua.
G= minreal(zpk([], [-1 -2], 5))
%step(G)
Tm= 0.09

%% FT discreta de lazo abierto Gd(s) del sistema con Z0H a la entrada
%% y el tiempo de muestreo asignado.
Gd= minreal(c2d(G, Tm, 'zoh'))

%% Mapa de polos y ceros del sistema continuo.
%figure(1);
%pzmap(G);

% Mapa de polos y ceros del sistema discreto.

%figure(2)
%zero(Gd)
%pole(Gd)
%pzmap(Gd)

% Qué ocurre con el mapa si se multiplica por 10 el periodo de muestreo?
%Gd1= c2d(G, 10*Tm, 'zoh');
%figure(3);
%pole(Gd1)
%zero(Gd1)
%pzmap(Gd1)
%%%%% CONCLUSIÓN
% Los polos se desplazan más hacia adentro del círculo unitario (hacia
% la izquierda).

%%%%% CONCLUSIÓN
%% Sistema inestable ya que los polos se encuentran dentro
%% del círculo unitario

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

Kp= dcgain(Gd);
Fd= feedback(Gd, 1)
step(Fd, 20)

% Almaceno los valores de la respuesta al escalón
[y, t]= step(Fd);
s_response= [y, t];

% Obtengo los valores de amplitud de la respuesta al escalón.
val_amplitud= s_response(:, 1);

% Realizo la diferencia de la ganancia
e= zeros(length(val_amplitud), 1);

for i= 1:1:length(val_amplitud)
   e(i)= Kp-val_amplitud(i);
end

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
Gd1= c2d(G, 10*Tm, 'zoh');

% Lugar de raíces (Se obtiene con realimentación unitaria).
%rlocus(G);
%title(sprintf("Lugar de raíces G"));
rlocus(Gd);
title(sprintf("Lugar de raíces Gd con Tm= %.2f", Tm));

%disp('Estabilidad relativa si se aumenta 10 veces el tiempo de muestreo')
%title(sprintf("Lugar de raíces Gd con Tm= %.2f", 10*Tm));
%rlocus(Gd1);
