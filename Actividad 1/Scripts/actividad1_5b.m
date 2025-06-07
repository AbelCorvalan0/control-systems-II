clear all; close all; clc;

data= readtable('../Resources/Curvas_Medidas_Motor_2025_v.xls');
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t).
t  = dataT(:,1);
wr = dataT(:,2);
ia = dataT(:,3);
v  = dataT(:,4);

figure(1); plot(t, wr); title('Angular velocity \omega_{R}'); 
ylabel('Angular velocity [rad/seg]'); xlabel('Time [seg]');
grid on; xlim([0, 0.7]);
figure(2); plot(t, ia); title('Current i_{a}(t)'); 
ylabel('Current [Ampere]'); xlabel('Time [seg]');
grid on;
figure(3); plot(t, v); title('Tension V'); 
ylabel('Voltage [V]'); xlabel('Time [seg]');
grid on; ylim([0, 3])