clc; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definición de parámetros físicos y condiciones iniciales

% Parámetros físicos del sistema
m = 0.1;         % Masa del péndulo (kg)
Fricc = 0.1;     % Coeficiente de fricción lineal (N·s/m)
long = 0.6;      % Longitud del péndulo (m)
g = 9.8;         % Aceleración gravitatoria (m/s²)
M = 0.5;         % Masa del carro (kg)

% Parámetros de simulación
h = 0.0001;              % Paso de integración (s)
tiempo = (50/h);         % Número total de iteraciones para 50 segundos

% Aceleraciones iniciales
p_pp = 0; tita_pp = 0;

% Condiciones iniciales: posición angular del péndulo
fi(1) = .1;       % Ángulo inicial del péndulo (rad)
color = 'r';      % Color para los gráficos

% Otros casos posibles comentados:
% fi(1)=.5; color='g';
% fi(1)=.8; color='b';
% fi(1)=-pi; color='+-b';
% fi(1)=pi; color='.-r';

omega(1) = 0;     % Velocidad angular inicial (rad/s)
p_p(1) = 0;       % Velocidad del carro inicial (m/s)
u(1) = 0;         % Control inicial
p(1) = 0;         % Posición inicial del carro
i = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linealización del sistema en el equilibrio inestable

% Versión linealizada del sistema alrededor del equilibrio inestable (phi = 0)
% Vector de estado: x = [p; p_p; fi; omega]

Mat_A = [ 0        1              0              0;
          0   -Fricc/M        -m*g/M            0;
          0        0              0              1;
          0  Fricc/(long*M)  g*(m+M)/(long*M)   0 ];

Mat_B = [ 0;
          1/M;
          0;
         -1/(long*M) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Verificación de controlabilidad

Mat_M = [Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B]; % Matriz de controlabilidad
rank(Mat_M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diseño del controlador LQR

% Definición de los pesos para el funcional de costo J = ∫(x'Qx + u'Ru)dt
Q = diag([.002 1.0 .10 .01]);  % Penaliza principalmente la velocidad del carro y el ángulo
R = 50;                        % Penaliza la magnitud del control

% Construcción del Hamiltoniano (matriz extendida de Riccati)
H = [ Mat_A       -Mat_B*inv(R)*Mat_B';
     -Q           -Mat_A' ];

[n_, va] = size(H);

% Cálculo de autovalores y autovectores del Hamiltoniano
[V,D] = eig(H); %V:autovectores, cada columna es un autovector, D:autovalores
MX1X2 = [];

% Se seleccionan las soluciones estables (autovalores con parte real negativa)
for(ii=1:n_)
    if real(D(ii,ii)) < 0
        MX1X2 = [MX1X2 V(:,ii)];
    end
end

% Separación en bloques
MX1 = MX1X2(1:n_/2,:); %M
MX2 = MX1X2(n_/2+1:end,:); %PM

% Cálculo de la solución P de Riccati (P simétrica definida positiva)
P = real(MX2*inv(MX1));

% Ganancia óptima de realimentación de estados
K = inv(R)*Mat_B'*P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulación del sistema no lineal con control LQR

estado = [p(1); p_p(1); fi(1); omega(1)];
J_(1) = 0; V_(1) = 0;

while(i < (tiempo + 1))
    estado = [p(i); p_p(i); fi(i); omega(i)];

    V_(i) = estado' * P * estado; % Función de Lyapunov

    fi1 = estado(3); % Ángulo actual
    
    % Corrección del ángulo para mantenerlo dentro de [-pi, pi]
    if abs(fi1) > pi
        while abs(fi1) > 2*pi
            fi1 = sign(fi1) * (abs(fi1) - 2*pi);
        end
        if fi1 < -pi
            fi1 = 2*pi + fi1;
        end
        if fi1 > pi
            fi1 = -2*pi + fi1;
        end
    end
    estado(3) = fi1;

    % Control óptimo (con saturación posterior)
    u(i) = -K * estado;

    % Funcional acumulado del costo
    J_(i+1) = J_(i) + (estado' * Q * estado + u(i)' * R * u(i)) * h;

    % Saturación del control
    u(i) = min(100, u(i));
    u(i) = max(-100, u(i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integración del sistema no lineal

    % Dinámica no lineal del sistema (Ecuaciones de Newton)
    p_pp = (1/(M + m)) * (u(i) - m * long * tita_pp * cos(fi(i)) + m * long * omega(i)^2 * sin(fi(i)) - Fricc * p_p(i));
    tita_pp = (1/long) * (g * sin(fi(i)) - p_pp * cos(fi(i)));

    % Integración por método de Euler
    p_p(i+1) = p_p(i) + h * p_pp;
    p(i+1) = p(i) + h * p_p(i);

    omega(i+1) = omega(i) + h * tita_pp;
    fi(i+1) = fi(i) + h * omega(i);

    estado_p = Mat_A * estado + Mat_B * u(i); % Evolución del modelo lineal (opcional)

    i = i + 1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gráficos

V_(i) = estado' * P * estado;
t = 0:tiempo; t = t * h;

% Figura 1: evolución temporal
figure(1); hold on;
subplot(3,2,1); plot(t, p, color); grid on; title('Posición carro');
subplot(3,2,3); plot(t, p_p, color); grid on; title('Velocidad carro');
subplot(3,2,2); plot(t, fi, color); grid on; title('Ángulo');
subplot(3,2,4); plot(t, omega, color); grid on; title('Velocidad angulo');
subplot(3,1,3); plot(t(1:end-1), u, color); grid on; title('Acción de control'); xlabel('Tiempo en Seg.');

% Figura 2: análisis cualitativo y funciones de energía
figure(2); hold on;
subplot(2,2,1); plot(fi, omega, color); grid on; xlabel('Ángulo'); ylabel('Velocidad angular');
subplot(2,2,2); plot(p, p_p, color); grid on; xlabel('Posición carro'); ylabel('Velocidad carro');
subplot(2,2,3); semilogy(t, J_, color); title('Funcional de costos J(x,u)'); xlabel('Tiempo Seg.');
subplot(2,2,4); semilogy(t, V_, color); title('Funcion Liapunov V(x)'); xlabel('Tiempo Seg.');

% Guardado de figuras
figure(1); print('-dpng','-r300','Fig_1');
figure(2); print('-dpng','-r300','Fig_2');