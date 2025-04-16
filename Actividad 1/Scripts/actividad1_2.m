%% Sistemas de Control II - FCEFyN (UNC) 
%
% En el archivo Curvas_Medidas_RLC_2025.xls (datos en la hoja 1 y etiquetas
% en la hoja 2) están las series de datos que sirven para deducir los 
% valores de R, L y C del circuito. Emplear el método de la respuesta al 
% escalón, tomando como salida la tensión en el capacitor.
%
%%

clear all; close all; 
clc;

%% Processing table
% Read table.
data= readtable('../Resources/Curvas_Medidas_RLC_2025.xls');
%t= csvread('../Resources/Curvas_Medidas_RLC_2025.xls',);
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t), capacitor voltage(Vcap) and current(I) array.
t= dataT(:,1);
I= dataT(:,2);
Vcap= dataT(:,3);

% Plot capacitor voltage and current

figure(1)

subplot(2, 1, 1);
plot(t, Vcap, 'b-', 'LineWidth', 1);
title('Capacitor voltage V_c');
ylabel('Voltage [V]');
xlabel('Time [seg]');
ylim([-15, 15]);
grid on;
hold on;

subplot(2, 1, 2);
plot(t, I, 'r-', 'LineWidth', 1);
title('Current I');
ylabel('Current [A]');
xlabel('Time [seg]');
grid on;
hold off;

% We've a step response with feature overdamped
% due to different real poles.

% Solve Equations 

% Define symbolic variables.
syms s;
syms R L C;
syms I Ve Vc;

% Define equations.
eq1= s*I== (1/L)*(Ve-Vc-R*I);
eq2= s*Vc== (1/C)*I;

% Obtain Vc, Ve solution.
s1= solve([eq1, eq2], [Vc, Ve]);
disp(s1.Vc);
disp(s1.Ve);

% Build transfer function.
TF= (s1.Vc)/(s1.Ve)
% We have a Second Order System.

% Apply The Second-Order Systems with differents poles method 

% y(t1)= K*(1+ ((T3-T1)/(T1-T2))*exp(-(t1/T1))- (((T3-T2)/(T1-T2))*exp(-t1/T2)));
% y(2*t1)= K*(1+ ((T3-T1)/(T1-T2))*exp(-(2*t1/T1))- (((T3-T2)/(T1-T2))*exp(-(2*t1)/T2)));
% y(3*t1)= K*(1+ ((T3-T1)/(T1-T2))*exp(-(3*t1/T1))- (((T3-T2)/(T1-T2))*exp(-(3*t1)/T2)));

% Obtain maximium value from t (RLC circuit measure).
max_t= max(t);

% Define step amplitude.
stepAmplitude= 1;

% Take time and voltage arrays
% Time array.
t0= t;
%Vc array.
y= dataT(:,3);

% Chen Method application
% Take three point to apply Chen Method.
% First point time.
%i= 1251;
i= 1074;
t_inic= t0(i);
% Define step.
h= 500;
% h= 55;
% Obtain y1.
t_t1= t0(i);
y_t1= y(i)
% Obtain y2.
t_2t1= t0(i+h);
y_2t1= y(i+h)
% Obtain y3.
t_3t1= t0(i+(2*h));
y_3t1= y(i+(2*h))

% Normalize gain
% Add abs(y(end)), because 'end' take value= -12V (Alternate input signal).
K= abs(y(end))/(stepAmplitude);

% Calculating k1, k2, k3.
k1= ((y_t1)/K)-1;
k2= ((y_2t1)/K)-1;
k3= ((y_3t1)/K)-1;

% Calculating b, alfa1, alfa2.
b= 4*(k1^(3))*k3- 3*(k1^2)*(k2^2)- 4*(k2^3)+ (k3^2)+ 6*k1*k2*k3;
alfa1= (k1*k2+ k3- sqrt(b))/(2*((k1^2)+ k2));
alfa2= (k1*k2+ k3+ sqrt(b))/(2*((k1^2)+ k2));

% Calculating Beta.
% beta= (2*(k1^3)+ 3*k1*k2+ k3- sqrt(b))/(sqrt(b));
% Alternative
beta= (k1+alfa2)/(alfa1-alfa2);

% Calculating the estimates of the time constants T1, T2 and T3.
% T1_e= -(t_t1)/(log(alfa1))
% T2_e= -(t_t1)/(log(alfa2)) %imaginary from alfa2
% deltaTime= 5.5e-4
T1_e= -(5.5e-4)/(log(alfa1));
T2_e= -(5.5e-4)/(log(alfa2));
T3_e= beta*(T1_e- T2_e)+ T1_e;

% Add estimates of the time constants T1, T2, T3.
ii= 1;
T1(ii)= T1_e;
T2(ii)= T2_e;
T3(ii)= T3_e;

% Enhancing estimation accuracy.
T1_e= sum(T1/ length(T1));
T2_e= sum(T2/ length(T2));
T3_e= sum(T3/ length(T3));

% Build Transfer Function
s= tf('s');
sys= (K)/((T1_e*s+1)*(T2_e*s+1));
sys1= (K*exp(-s*0.01))/((T1_e*s+1)*(T2_e*s+1))
% [num, den]= tfdata(sys, 'v');

% Plot approximated step response
% figure(2);
% %%%%% FIX
% Obtain step response with delay.
% [ys, ts]= step(sys*exp(-s*0.01), 'r-', 0.16);
% [ys, ts]= step(sys*exp(-s*0.01), 'r-');
% ylim([0, 13]);
% title('Approximation');
% grid on;
% figure(3);
% plot(ts, ys, t, Vcap);
% plot(t, Vcap);
% plot(ts, ys);
% title('Approximation vs Measurement');
% xlabel('Amplitude[V]')
% ylabel('Time[seg]')
% grid on;
% ylim([0, 13]);
% xlim([0.008, 0.014]);
% legend('Approximated', 'Measured');
% hold off;
% 
% size(ys)
% size(ts)
% 
% tab= [real(ts), real(ys)];
% Generate .csv to approximated step response 
% filename= 'Generated Data/approximation.csv';
% headers= {'Time [Seconds]', 'Amplitude [V]'};
% fi= fopen(filename, 'w');
% fprintf(fi, '%s,%s\n', headers{:});
% fclose(fi);
% writematrix(tab, filename, 'WriteMode', 'append');
% 
% Obtain Parameters R, L, C 
% 
% Assume resistance value.
% R= 220;
% 
% LC= den(1, 1);
% RC= den(1, 2);
% 
% C= real(RC/R)
% L= real(LC/C)
% 
% Build State-Space 
% 
% A= [-R/L -1/L; 1/C 0];
% B= [1/L; 0];
% c= [R 0];
% D= 0;
% 
% Function state space model
% sys= ss(A, B, c, D);
% 
% Input signal 
% 
% Parameters
% h= 1e-5;
% t= 0: h: 29.999e-3; 
% u= zeros(size(t));
% var= 0;
% tinit= 5e-3;
% T= 10e-3;
% 
% for i=1:length(t)
%     if t(i)<= tinit
%         u(i)= 0;
%     elseif t(i)> tinit && var<= (T/h)
%         u(i)= 12;
%         var= var+1;
%     elseif t(i)> tinit && var> (T/h)
%         u(i)= -12;
%         var= var+1;
%     end
% 
%     if var> (2*(T/h))
%         var= 0;
%     end
% end
% 
% Plot input signal 
% figure(4);
% plot(t, u);
% grid on;
% xlim([0, 0.030]);
% ylim([-13, 13]);
% xlabel('Time [Seconds]');
% ylabel('Voltage [V]');
% title('Input Signal');
% 
% State-Space Model Simulation
% figure(5);
% lsim(sys, u, t);
% d= lsim(sys, u, t);
% size(z)
% size(t')
% data= [t', d];
% grid on;
% 
% %%%% Generate .csv file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% csvwrite('simulation1.csv', data);
% filename= 'Generated Data/simulation1.csv';
% 
% headers= {'Time [Seconds]', 'Amplitude [V]'};
% fi= fopen(filename, 'w');
% fprintf(fi, '%s,%s\n', headers{:});
% fclose(fi);
% writematrix(data, filename, 'WriteMode', 'append');