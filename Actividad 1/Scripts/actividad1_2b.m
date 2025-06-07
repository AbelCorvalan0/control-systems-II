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
Is= dataT(:,2);
Vcap= dataT(:,3);

% Plot capacitor voltage and current

figure(1)
subplot(2, 1, 1);
plot(t, Vcap, 'b-', 'LineWidth', 1);
title('Capacitor Voltage v_{c}(t)');
ylabel('Voltage [V]');
xlabel('Time [seconds]');
ylim([-15, 15]);
grid on;
hold on;

subplot(2, 1, 2);
plot(t, Is, 'r-', 'LineWidth', 1);
title('Current i_{a}(t)');
ylabel('Current [A]');
xlabel('Time [seconds]');
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

% y(t1)   = K*(1+ ((T3-T1)/(T1-T2))*exp(-(t1/T1))- (((T3-T2)/(T1-T2))*exp(-t1/T2)));
% y(2*t1) = K*(1+ ((T3-T1)/(T1-T2))*exp(-(2*t1/T1))- (((T3-T2)/(T1-T2))*exp(-(2*t1)/T2)));
% y(3*t1) = K*(1+ ((T3-T1)/(T1-T2))*exp(-(3*t1/T1))- (((T3-T2)/(T1-T2))*exp(-(3*t1)/T2)));

% Obtain maximium value from t (RLC circuit measure).
max_t= max(t);

% Define step amplitude.
stepAmplitude= 12;
K= abs(Vcap(end));

% Take time and voltage arrays
% Time array.
t0 = t;
%Vc array.
y  = dataT(1001:end,3);
t0 = t0(1001:end,1);
% Vcap1 = Vcap(1001:end, 1);

% Chen Method application
% Take three point to apply Chen Method.
% First point time.
i= 35;
t_inic= t0(i)
disp('Disp')
% % Define step.
% h= 200;
% Obtain y1.
t_t1= t0(i)
y_t1= y(i)
% Obtain y2.
t_2t1= t0(2*i)
y_2t1= y(2*i)
% Obtain y3.
t_3t1= t0(3*i)
y_3t1= y(3*i)

delayTime= 0.01;

%%
% % Normalize gain
% % Add abs(y(end)), because 'end' take value= -12V (Alternate input signal).
% K= abs(y(end))/(stepAmplitude);
% 
% % Calculating k1, k2, k3.
% k1= ((y_t1)/K)-1;
% k2= ((y_2t1)/K)-1;
% k3= ((y_3t1)/K)-1;
% 
% % Calculating b, alfa1, alfa2.
% b= 4*(k1^(3))*k3- 3*(k1^2)*(k2^2)- 4*(k2^3)+ (k3^2)+ 6*k1*k2*k3;
% alfa1= (k1*k2+ k3- sqrt(b))/(2*((k1^2)+ k2));
% alfa2= (k1*k2+ k3+ sqrt(b))/(2*((k1^2)+ k2));
% 
% % Calculating Beta.
% % beta= (2*(k1^3)+ 3*k1*k2+ k3- sqrt(b))/(sqrt(b));
% % Alternative
% beta= (k1+alfa2)/(alfa1-alfa2);
% 
% % Calculating the estimates of the time constants T1, T2 and T3.
% % T1_e= -(t_t1)/(log(alfa1))
% % T2_e= -(t_t1)/(log(alfa2)) %imaginary from alfa2
% % deltaTime= 5.5e-4
% T1_e= -(5.5e-4)/(log(alfa1));
% T2_e= -(5.5e-4)/(log(alfa2));
% T3_e= beta*(T1_e- T2_e)+ T1_e;
% 
% % Add estimates of the time constants T1, T2, T3.
% ii= 1;
% disp('Constantes de tiempo: ')
% T1(ii)= T1_e;
% T2(ii)= T2_e;
% T3(ii)= T3_e;
% 
% % Enhancing estimation accuracy.
% T1_e= sum(T1/length(T1))
% T2_e= sum(T2/length(T2))
% T3_e= sum(T3/length(T3))
% 
% % Build Transfer Function
% s= tf('s');
% sys= (K)/((T1_e*s+1)*(T2_e*s+1));
% % sys= (K*(T3_e*s+1))/((T1_e*s+1)*(T2_e*s+1));
% sys1= (sys*exp(-s*0.01))
% % [num, den]= tfdata(sys, 'v');
%%

sys= chenMethod(y_t1, y_2t1, y_3t1, K, stepAmplitude, delayTime, t_t1, 2)

% Plot approximated step response
% Obtain step response with delay.
figure(2)
[ys, ts]= lsim(sys, dataT(:, 4), t);
plot(t, Vcap, 'b-'); hold on;
%plot(ts, ys, 'r-');
xlim([0, 0.025]); ylim([0, 13]);
xlabel('Time [seconds]'); ylabel('Voltage [V]');
plot(t_t1, y_t1, 'o'); hold on;
plot(t_2t1, y_2t1, 'o'); hold on;
plot(t_3t1, y_3t1, 'o'); hold on;
% legend('Measured', 'Approximated');
legend('Measured');
% title('Approximation vs Measurement');
title('Step Response');
grid on;
tab= [ts, ys];

%% Generate .csv to approximated step response 
filename= 'Generated Data/approximation.csv';
headers= {'Time [Seconds]', 'Amplitude [V]'};
fi= fopen(filename, 'w');
fprintf(fi, '%s,%s\n', headers{:});
fclose(fi);
writematrix(tab, filename, 'WriteMode', 'append');

%% Obtain Parameters R, L, C.
% Capacitor value calculate.
disp('Voltage')
Vcap(1145)
Vcap(1208)
disp('Time')
t(1145)
t(1208)
disp('Current')
Is(1145)

deltaV= (Vcap(1145)-Vcap(1208))/(t(1145)-t(1208))
it= Is(1145);

C= it/(deltaV)
% FT= (1/(LC))/(s^2+s*(R/L)+(C/L))
sys
[num, den]= tfdata(sys, 'v');
% FTn= minreal(sys);
% [num, den]= tfdata(FTn, 'v');

L= den(1)/C
R= den(2)/C