close all; clear all;
clc;
%% Tarea 3 - Laboret.



%% 1
%% Defino parámetros del sistema.
m= 2;

mnom=m % masa nominal
m1=0.9*mnom % correr de nuevo el código de simulación, dibujo y analisis
m2=1.1*mnom % ídem

b= 0.3;
l= 1;
G= 10;
delta= 180 % [grados].
p_triple= -4;

% linmod()
% computes the linear state-space model of the system of ordinary 
% differential equations represented in the model mdl by linearizing
% each block one by one. Inport and Outport blocks in the model 
% represent the system inputs and outputs.
[A, B, C, D]= linmod('pendulo_mod_tarea', delta*(pi/180))
disp('Los autovalores del sistema son: ')
eig(A)
disp('El rango de la matriz controlabilidad es: ')
rank(ctrb(A, B))

%% Torque de equilibrio.

uf= m*G*l*sin(delta*pi/180)
uf1= m1*G*l*sin(delta*pi/180)
uf2= m2*G*l*sin(delta*pi/180)

%% 2
%% Introducir integral del error.
% Ampliación de sistema.
disp('Sistema ampliado: ')
Aa= [ [A; C] zeros(3, 1)];

Ba= [ B ;
      0 ];

% Autovalores del sistema ampliado.
disp('Autovalores de la matriz ampliada: ')
eig(Aa)

% Verificar rango de la matriz controlabilidad.
disp('Rango de matriz controlabilidad')
rank(ctrb(Aa, Ba));
% verficación de controlabilidad del sistema.

%% 3
%% Diseñar por asignación de polos un controlador 
% Implementar acker().

% Se debe asignar un polo triple
p= -4;
K= acker(Aa, Ba, [p p p])
% Valores individuales de ganacia.
k1= K(1);
k2= K(2);
k3= K(3);

% Polos de lazo cerrado
disp('Polos de lazo cerrado: ')
eig(Aa-Ba*K)

% Tiempo de respuesta calculado
disp('Tiempo de respuesta calculado: ')
tscalc= 7.5/(-p)

%% 4
%% Simulación péndulo con PID.
% punto inicial origen
% velocidad nula y referencia delta= 180°.

sim('pendulo_pid_tarea')

% Plots
% Dibujar: 1- Salida.
%          2- Plano de fases.
%          3- Torque total.
%          4- Acción integral.

figure(1); plot(tout, yout);
grid on; title('Salida');

figure(2), plot(yout,velocidad) %plano de fase
grid on, title('Plano de fases');

figure(3), plot(tout,torque) % torque total
grid on, title('Torque');

figure(4), plot(tout,-accint) % acción integral
grid on, title('Accion integral')

%% Cálculos a comparar con plots.
disp('Cálculos: ')
ymax=max(yout) % máximo valor de salida
S=(ymax-delta)/delta*100 % sobrepaso en (%)
erel=(delta-yout)/delta; %error relativo
efinal=erel(end) % error final, debe ser cero
ind=find(abs(erel)>.02); % índice elementos con error relativo absoluto menor a 2%
tss=tout(ind(end)) % tiempo de establecimiento (ultimo valor del vector)
yte=yout(ind(end)) % salida al tiempo ts
uf=torque(end) % torque final
Intf=-accint(end) % acción integral final


% verificar sobrepaso. 
% Tiempo de establecimiento real vs calculado.

m_nom = 2;
variaciones_masa = [1.8, 2, 2.2];
colores = ['r', 'g', 'b']; % colores para las 3 curvas
labels = {'1.8 m', '2 m', '2.2 m'};

figure(10); hold on; title('Salida \theta(t)')
figure(11); hold on; title('Plano de fases ')
figure(12); hold on; title('Torque aplicado')
figure(13); hold on; title('Acción integral')

for i = 1:3
    m = variaciones_masa(i) * m_nom;
    sim('pendulo_pid_tarea');
    
    figure(10)
    plot(tout, yout, 'DisplayName', labels{i})
    
    figure(11)
    plot(yout, velocidad, 'DisplayName', labels{i})
    
    figure(12)
    plot(tout, torque, 'DisplayName', labels{i})
    
    figure(13)
    plot(tout, -accint, 'DisplayName', labels{i})
end

for f = 10:13
    figure(f)
    legend
    xlabel('Tiempo [s]')
    grid on
end

%%% Robustez

m_nom = 1;
variaciones_masa = [0.9, 1.0, 1.1];
resultados = zeros(3, 4); % Columnas: sobrepaso, tss, uf, Intf

for i = 1:3
    m = variaciones_masa(i) * m_nom;
    sim('pendulo_pid_tarea')

    ymax = max(yout);
    S = (ymax - delta)/delta * 100;
    erel = (delta - yout)/delta;
    ind = find(abs(erel) > 0.02);
    tss = tout(ind(end));
    uf = torque(end);
    Intf = -accint(end);

    resultados(i,:) = [S, tss, uf, Intf];
end

% Mostrar tabla
fprintf('\nTabla de robustez (masa variable):\n');
fprintf('Masa\t\tSobrepaso (%%)\tTs (s)\t\tUf\t\tIntf\n');
for i = 1:3
    fprintf('%.2f\t\t%.2f\t\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
        variaciones_masa(i)*m_nom, resultados(i,1), resultados(i,2), resultados(i,3), resultados(i,4));
end
%%
m_nom = 2;
variaciones_masa = [1.8, 2, 2.2];
colores = ['r', 'g', 'b']; % colores para las 3 curvas
labels = {'1.8 m', '2 m', '2.2 m'};

figure(10); hold on; title('Respuesta angular \theta(t)')
figure(11); hold on; title('Plano de fases (\theta, \dot{\theta})')
figure(12); hold on; title('Torque aplicado')
figure(13); hold on; title('Acción integral acumulada')

for i = 1:3
    m = variaciones_masa(i) * m_nom;
    sim('pendulo_pid_tarea');
    
    figure(10)
    plot(tout, yout, 'DisplayName', labels{i})
    
    figure(11)
    plot(yout, velocidad, 'DisplayName', labels{i})
    
    figure(12)
    plot(tout, torque, 'DisplayName', labels{i})
    
    figure(13)
    plot(tout, -accint, 'DisplayName', labels{i})
end

for f = 10:13
    figure(f)
    legend
    xlabel('Tiempo [s]')
    grid on
end
