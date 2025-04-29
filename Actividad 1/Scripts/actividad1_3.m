clear all; close all; clc;
%%% Obtain Parameters R, L, C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R= 119.6863;
C= 3.9434e-06;
L= 0.0011;

%% Simulate State-Space 
num= [C 0];
den= [L*C C*R 1];
s= tf('s');
G= tf(num, den)
%[A, B, c, D]= tf2ss(num, den);
p= pole(G);
% Fastest dynamic pole.
p(1)
% Slower dynamic pole.
p(2)

vin= 12;
tint = log(0.95)/p(1);
tsim = log(0.05)/p(2);
T= 10e-3;
tsim= 20*T;
h= tsim/tint;
t= linspace(0, tsim, h);
u= linspace(0, 0, h);
var=0;

% Hacer variable las operaciones que se repiten
% por ejemplo round(T/tint)
h10m= round(4*T/tint);

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

%% Operation point.
% Initial values x(0)
Il(1)  = 0;        % Current.
Vrl(1) = 0;        % Voltage through R. ??
y(1)   = 0;        % Output.

Xop    = [0, 0]';
x      = [Il(1) Vrl(1)]'; 

%% State Space.
A= [-R/L -1/L; 1/C 0];   % State matrix.
B= [1/L; 0];             % Input matrix.
c= [R 0];                % Output matrix.
D= 0;                    % Direct transmission matrix.

for i=1: h-1
    if t(i)<= 0.01
        u(i)= 0;
    elseif t(i)> T && var<=(h10m)
        u(i)= 12;
        var= var+1;
    elseif t(i)>= T && var>(h10m)
        u(i)= -12;
        var= var+1;
    end
    if var>(2*h10m)
        var=0;
        u(i)= 0;
    end
    xp= A*(x-Xop)+B*u(i);
    x = x+ xp*tint;
    Y= c*(x-Xop)+ D*u(i);
 
    ii= i+1;
    y(ii)   = Y(1);
    Il(ii)  = x(1);
    Vrl(ii) = x(2);

end

data= readtable('../Resources/Curvas_Medidas_RLC_2025.xls');
%t= csvread('../Resources/Curvas_Medidas_RLC_2025.xls',);
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t), capacitor voltage(Vcap) and current(I) array.
t1= dataT(:,1);
Is= dataT(:,2);
Vcap= dataT(:,3);
%max(Is)

figure(1)
subplot(3, 1, 1);
plot(t, u); grid on;
xlim([0, 0.08]); ylim([-13, 13]);
xlabel('Time [seconds]'); ylabel(['Voltage [V]']);
title('Input Signal u(t)');
subplot(3, 1, 2);
plot(t, Il); grid on; hold on;
xlim([0, 0.08]); ylim([-0.2, 0.2]);
xlabel('Time [seconds]'); ylabel('Current [mA]')
title('Approximated Current Signal I_{s}(t)')
subplot(3, 1, 3)
plot(t1, Is); grid on;
xlim([0, 0.08]); ylim([-0.2, 0.2]);
xlabel('Time [seconds]'); ylabel(['Current [mA]'])
title('Measured Current Signal I_{s}(t)')

%%%%% Generate .csv file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%csvwrite('simulation1.csv', data);
% filename= 'Data Generated/simulation1.csv';
% 
% headers= {'Time [Seconds]', 'Amplitude [I]'};
% fi= fopen(filename, 'w');
% fprintf(fi, '%s,%s\n', headers{:});
% fclose(fi);
% writematrix(data, filename, 'WriteMode', 'append');