%% Sistemas de Control II - FCEFyN (UNC) 
%
% Asignar valores a R= 220Ohm, L= 500mHy, y C= 2,2uF. Obtener simulaciones 
% que permitan estudiar la dinámica del sistema, con una entrada de tensión 
% escalón de 12V, que cada 10ms cambia de signo.
%
%%

clc; clear all;
close all;

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

%% Design Values

R= 220;
L= 500e-3;
C= 2.2e-6;

%% Operation point.
Il(1)  = 0;
Vcl(1) = 0; 
y(1)   = 0;

Xop    = [0 0]';
x      = [Il(1) Vcl(1)]';

%% State variables.
A= [-R/L -1/L; 1/C 0];   % State matrix
B= [1/L; 0];             % Input matrix
c= [R 0];                % Output matrix
D= 0;                    % Direct transmission matrix

%% Function state space model
sys= ss(A, B, c, D);

[num, den]= ss2tf(A, B, c, D);
% Build transfer function.
G= tf(num, den)

%% System Analysis
% Obtain poles.
poles= pole(G)
damp(G)
[wn, zeta]= damp(G);
wn(1);
zeta(1);
% Analysis of zeta
% zeta= 0.2310 
% 0<zeta<1 underdamped temporal response.

% Calculate damped angular frequency (wd).
wd= wn(1)*sqrt(1-(zeta(1)^2))
%wd= 927.7343 rad/seg.

% Calculate natural period Tn.
Tn= (2*pi)/wd;

% Integration time. 
tint= Tn/100

% Simulation time.
% tsim= 3*log(0.05)/(real(poles(1)))
tsim= 10000;     % Arbitrarily determinated.
% Three times of time which exp(lambda2*t) reachs 5%.

h= tsim/tint;

%% Input signal
% Step time.
h= 0.0001;  % Step of 0.1ms.

t= 0: h: 0.1999;    % Time vector ranged 0s until 200ms (step= 0.1ms). 
                    % (defined by h). We have 2000 positions.
u= zeros(size(t));  % Input signal.
variable= 0; % Counts 10ms/h= 100 steps.
             % Switch voltage each 10ms.

for i = 1:length(t)
    
    if t(i) <= 0.02
        u(i)= 0; % Firs 20ms takes zero (0) value.
    elseif t(i) > 0.02 && variable<=(0.01/h)  % Take 100 steps= 10ms.
        u(i) = 12;
        variable=variable+1;
    elseif t(i) > 0.02 && variable>(0.01/h)
        u(i) = -12;
        variable=variable+1;
    end
    
    if variable>(2*(0.01/h)) % Reset
        variable=0;
    end
end

%% Plot input signal (u)

figure(1);
plot(t, u);
grid on;
%xlim([0, 0.02]);
ylim([-13, 13]);
xlabel('Time [Seconds]');
ylabel('Voltage [V]');
title('Input Signal');


%% State-Space Model Simulation
figure(2);
lsim(sys, u, t);
xlim([0, 0.06])
z= lsim(sys, u, t);
%size(z)
%size(t')
data= [t', z];
grid on;

%% Generate .csv file
%csvwrite('simulation1.csv', data);
filename= 'Generated Data/simulation1.csv';

headers= {'Time [Seconds]', 'Amplitude [V]'};
fi= fopen(filename, 'w');
fprintf(fi, '%s,%s\n', headers{:});
fclose(fi);
writematrix(data, filename, 'WriteMode', 'append');