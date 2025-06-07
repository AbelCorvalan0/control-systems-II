clear all; close all; clc;

%% Operation point.
% % Initial values x(0)
Ia(1)    = 0;        % 
wr(1)    = 0;        % 
theta(1) = 0;        %

y(1)     = 0;        % Output.

Xop    = [0, 0, 0]';
x      = [Ia(1) wr(1) theta(1)]'; 

Laa= 366e-6;    % [Hy].
J  = 5e-9;      % [kg.m^2].
Ra = 55.6;      % [Ohm].
Bm = 0;         % [N.m/(rad/seg)].
Ki = 6.49e-3;   % 
Km = 6.53e-3;   % 

% Laa  = 0.002511;
% Ki  = 3.765;
% J   = 0.02056;
% Km  = 0.2485;
% Bm  = 0.026;
% Ra  = 2.415;

A= [-Ra/Laa  -Km/Laa  0;
    Ki/J     -Bm/J    0;
     0         1      0];

B= [1/Laa   0 ;
      0   -1/J;
      0     0 ];

C= [0 1 0];

D= [0 0];

%% Simulation
tint= 10e-7;

%% Simulation time calculate.

ts= 5;
h= ts/tint;

t  = 0: tint: ts;
%Va = linspace(0, 0, tint);
Va = zeros(1, length(t));
Tl = zeros(1, length(t));
Tlval= 0;
%Tlval= 1.162;
%Tlval= 0.0012;

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

for i=1: length(t)-1
    if t(i)>= 0
        Va(i)= 12;
    else 
        Va(i)= 0;
    end

    if t(i)>= 0.701 && t(i)<= 1.001 
        Tl(i)= Tlval;
    else
        Tl(i)= 0;
    end
   
    xp = A*(x-Xop) + B*[Va(i) Tl(i)]';
    x =  x + xp*tint;
    Y  = C*(x-Xop) + D*[Va(i) Tl(i)]';
    
    ii        = i+1;
    wr(ii)    = Y(1);
    Ia(ii)    = x(1);
    theta(ii) = x(3);
end
disp('W_{R \, max}= ')
wr(end)
disp('I_{a \, max}= ')
max(Ia)
 
disp(['Maximium torque calculating: T= Km*ia'])
Torque= Km*abs(max(Ia))

figure(1);
subplot(5, 1, 1);
plot(t, Va); title('Input Signal'); ylabel('Voltage [V]'); xlabel('Time [seg]');
xlim([0, 2.5]);
subplot(5, 1, 5);
plot(t, wr); title('Angular Velocity \omega_{R}'); ylabel('Angular Velocity [rad/seg^{2}]'); xlabel('Time [seg]');
xlim([0, 2.5]); grid on;
% subplot(5, 1, 3);
% plot(t, theta); title('Angular Position'); ylabel('Angular Posit [rad/seg]'); xlabel('Time [seg]');
% xlim([0, 2.5]);
subplot(5, 1, 2);
plot(t, Ia); title('Current i_{a}'); ylabel('Current [A]'); xlabel('Time [sec]');
xlim([0, 2.5]); grid on; %ylim([])
subplot(5, 1, 3);
plot(t, Tl); title('Torque Load T_{L}'); ylabel('Torque [Nm]'); xlabel('Time [sec]');
xlim([0, 2.5]); grid on; %ylim([])

subplot(5, 1, 4)
plot(t, theta); title('Angular position \theta_{t}'); 
ylabel('Angular position [rad]'); xlabel('Time [sec]')