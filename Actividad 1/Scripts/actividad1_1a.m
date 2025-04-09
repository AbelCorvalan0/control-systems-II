%% Sistemas de Control II - FCEFyN (UNC)
%
% Asignar valores a R= 220Ohm, L= 500mHy, y C= 2,2uF. Obtener simulaciones 
% que permitan estudiar la dinámica del sistema, con una entrada de tensión 
% escalón de 12V, que cada 10ms cambia de signo.
%%

clc; clear all;
close all;

%%
% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]
%%

R= 220;
L= 500e-3;
C= 2.2e-6;

% Operation point.
Il(1)  = 0;
Vcl(1) = 0;
y(1)   = 0;

Xop    = [0 0]';
x      = [Il(1) Vcl(1)]';

% State variables.
A= [-R/L -1/L; 1/C 0];   % State matrix
B= [1/L; 0];             % Input matrix
c= [R 0];                % Output matrix
D= 0;                    % Direct transmission matrix

% Function state space model
sys= ss(A, B, c, D);

[num, den]= ss2tf(A, B, c, D);
G= tf(num, den)

% Obtain poles.
poles= pole(G)
damp(G)
[wn, zeta]= damp(G);
wn(1);
zeta(1);
% Analysis of zeta
% zeta= 0.2310 
% 0<zeta<1 underdamped temporal response.

% Calculate damped angular frequency (\omega_{d}).
wd= wn(1)*sqrt(1-(zeta(1)^2))
%wd= 927.7343 rad/seg.

% Calculate natural period Tn.
Tn= (2*pi)/wd;

% Integration time. 
tint= 10*Tn

% Simulation time.
%tsim= 3*log(0.05)/(real(poles(1)))
tsim= 10000;
% Three times of time which exp(lambda2*t) reachs 5%.

h= tsim/tint;

%%%%%% Input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tiempo de integración
h= 0.0001;  % 0.1ms de paso

t= 0: h: 0.1999;    % Vector de tiempo de 0 a 200ms con paso de 0.1ms 
                   % (definido por h). Se tienen 2000 posiciones.
u= zeros(size(t));
variable= 0; % Cuenta 10ms/h= 100 pasos
             % Se cambia de estado cada 10ms

for i = 1:length(t)
    
    if t(i) <= 0.02
        u(i)= 0; %primeros 20ms vale cero
    elseif t(i) > 0.02 && variable<=(0.01/h)  % Realiza 100 pasos= 10ms.
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

figure(1);
plot(t, u);
grid on;
%xlim([0, 0.02]);
ylim([-13, 13]);
xlabel('Time [Seconds]');
ylabel('Voltage [V]');
title('Input Signal');


%%%% State-Space Model Simulation %%%%%%%%%%%%%%%%%
figure(2);
lsim(sys, u, t);
xlim([0, 0.06])
z= lsim(sys, u, t);
%size(z)
%size(t')
data= [t', z];
grid on;

%%%%% Generate .csv file %%%%%%%%%%%%%%%%%%%%
%csvwrite('simulation1.csv', data);
filename= 'Generated Data/simulation1.csv';

headers= {'Time [Seconds]', 'Amplitude [V]'};
fi= fopen(filename, 'w');
fprintf(fi, '%s,%s\n', headers{:});
fclose(fi);
writematrix(data, filename, 'WriteMode', 'append');