clear all; 
close all; clc;

%% Import data table.
% Read table.
data= readtable('../Resources/Curvas_Medidas_Motor_2025_v.xlsx');
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t), capacitor voltage(Vcap) and current(I) array.
t  = dataT(:, 1);
wr = dataT(:, 2);
I  = dataT(:, 3);
V  = dataT(:, 4);
tl = dataT(:, 5);
disp('Lenght of tl vector')
length(t)
disp('Maximium value time.')
max(t)

t(2)-t(1)

subplot(4, 1, 1);
plot(t, tl, 'LineWidth', 1.5)
xlabel('Tiempo [seg]')
ylabel('Torque [Nm]')
title('Perturbación (Torque) \T_{L}')
grid

%% Define system parameters.
%
Ra  = 2.258930051299405;
La  = 0.005026901184834716;
Km  = 0.2500481104997174;
Ki  = 0.25965987053759737;
Jm  = 0.0028472626983113334;
Bm  = 0.0014165170369840668;


%% Define State Space Model.

% State Matrix.
A= [-Ra/La -Km/La  0;
     Ki/Jm -Bm/Jm  0;
       0      1    0];

% Input matrix.
B= [1/La    0 ;
     0   -1/Jm;
     0      0 ];

% Output Matrix.
C= [0 0 1];

D= [0 0];

%% Initial conditions.
Xop      = [0 0 0]';
ia(1)    = 0;
wr(1)    = 0;
theta(1) = 0;

sys= ss(A, B, C, D)

%% Controllability verification.
% Controllability matrix calculation.
Co= ctrb(sys)
% It allows us if the system can be carried brought into a any state
% in a limited time.

% Range of controllability matrix.
%rank_Co= rank(Co);

if rank(Co) == size(A, 1)
    disp('Sistema controlable');
else
    error('Sistema no controlable');
end   

% Eigen values from A.
% eigenV= eig(A)  % display eigenvalues.
% We've two negative real part eigenvalues
% stable system.
% Besides we've third eigenvalue which is equal to zero (0).
% which makes the system marginally stable.

%% Expanded matrix.
Aa= [A  zeros(3,1);
    -C  0]

Ba= [B(:,1); 0]

%% 
%tl = ((1.15e-3)/2)*square(2*pi*(1/50)*t)+((1.15e-3)/2);

%% LQR Controller Design
% Define weighting matrices Q and R
Q = diag([0.5, 0.5, 100, 1000]);  % Weight on states (position error heavily weighted)
R = 100;      % Weight on control inputs

% Compute LQR gain matrix K.
[K, S, e] = lqr(Aa, Ba, Q, R);

poles= eig(Aa-Ba(:,1)*K)

%% Simulation time calculation.
tsim = 50;

%% Integration time.
tint = (log(.5)/(min(poles)))/4
% importante

simTime= 50;
h= 1e-4;
t= 0:h:(simTime-h);

% Build reference (pi/2, -pi/2).
reference = (pi/2)*square(2*pi*(1/10)*t);

%figure(1)
subplot(4, 1, 2);
plot(t,reference, 'b--','LineWidth', 1.5)
xlabel('Tiempo [seg]')
ylabel('Ángulo [rad]')
title('Señal de referencia \theta_i')
grid

tl = ((1.15e-3)/2)*square(2*pi*(1/50)*t)+((1.15e-3)/2);
% subplot(4, 1, 1);
% plot(t, tl, 'LineWidth', 1.5)
% xlabel('Tiempo [seg]')
% ylabel('Torque [Nm]')
% title('Perturbación (Torque) \T_{L}')
% grid

stateVector = [ia(1) wr(1) theta(1)]';
x= [ia(1) wr(1) theta(1)];
Xop = [0 0 0]';
zeta(1)  = 0;
integ(1) = zeta(1);
K

% torque begin at 0.701seg.
% 0.12 Nm of amplitude.

for i = 1:(simTime/h)

    zetaP = reference(i)-C*stateVector;
    zeta(i) = integ+zetaP*h;
    u(i)  = -K(1:3)*stateVector-K(4)*zeta(i);
    ia(i)    = x(1);
    wr(i)    = x(2);
    theta(i) = x(3);
    x1P = -Ra*x(1)/La-Km*x(2)/La+u(i)/La;
    x2P = Ki*x(1)/Jm-Bm*x(2)/Jm-tl(i)/Jm;
    x3P = x(2);
    xP = [x1P x2P x3P]';
    x = x+h*xP;
    stateVector = [ia(i) wr(i) theta(i)]';
    integ = zeta(i);
end

%figure(2)
subplot(4, 1, 2); hold on;
plot(t, theta, 'LineWidth', 1.5)
xlabel('Tiempo [seg]')
ylabel('Ángulo [rad]')
title('Salida \theta_r')
grid

%% 
% Junto con la salida del ángulo del motor, también se pueden visualizar la 
% corriente de armadura y la acción de control.
% figure(3)
subplot(4, 1, 3);
plot(t, ia, 'LineWidth', 1.5)
xlabel('Tiempo [seg]'); ylabel('Corriente [A]');
ylim([-0.25, 0.25])
title('Salida i_a')
grid

%figure(4)
subplot(4, 1, 4);
plot(t, u, 'LineWidth', 1.5)
xlabel('Tiempo [seg]'); ylim([-1.5, 1.5])
title('Acción de control $u(t)$')
grid