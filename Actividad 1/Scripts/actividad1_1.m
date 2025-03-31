clc; clear all;
close all;

R= 220;
L= 500e-3;
C= 2.2e-6;
    
A= [-R/L -1/L; 1/C 0];
B= [1/L; 0];
c= [R 0];
D= 0;

% Function state space model
sys= ss(A, B, c, D);

%%%%%% Input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters %%%
% Steps.
h= 1e-5;
% Time vector
t= 0: h: 2.499999e-2;
% Amplitude vector.
u= zeros(size(t));
var= 0;
% Init time.
tinit= 5e-3;
T= 1e-3;

for i=1:length(t)
    if t(i)<= tinit
        u(i)= 0;
    elseif t(i)> tinit && var<= (T/h)
        u(i)= 12;
        var= var+1;
    elseif t(i)> tinit && var> (T/h)
        u(i)= -12;
        var= var+1;
    end

    if var> (2*(T/h))
        var= 0;
    end
end

figure(1);
plot(t, u);
grid on;
xlim([0, 0.02]);
ylim([-13, 13]);
xlabel('Time [Seconds]');
ylabel('Voltage [V]');
title('Input Signal');


%%%% State-Space Model Simulation %%%%%%%%%%%%%%%%%
figure(2);
lsim(sys, u, t);
z= lsim(sys, u, t);
%size(z)
%size(t')
data= [t', z];
grid on;

%%%%% Generate .csv file %%%%%%%%%%%%%%%%%%%%
%csvwrite('simulation1.csv', data);
filename= 'Data Generated/simulation1.csv';


headers= {'Time [Seconds]', 'Amplitude [V]'};
fi= fopen(filename, 'w');
fprintf(fi, '%s,%s\n', headers{:});
fclose(fi);
writematrix(data, filename, 'WriteMode', 'append');