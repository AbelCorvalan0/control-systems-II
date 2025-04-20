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

s   = tf('s');
num = [Ki/(J*Laa)];
den = [1,((Bm*Laa + J*Ra))/J, (Laa*(Bm*Ra + Ki*Km))/J];

TF= tf(num, den)

% [A, B, C, D]= tf2ss(num, den);

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

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

for i=1: length(t)-1
    if t(i)>= 0.5
        Va(i)= 12;
    else 
        Va(i)= 0;
    end
    
    xp = A*(x-Xop) + B*[Va(i) Tl(i)]';
    x =  x + xp*tint;
    Y  = C*(x-Xop) + D*[Va(i) Tl(i)]';
    
    ii        = i+1;
    wr(ii)    = Y(1);
    Ia(ii)    = x(1);
    theta(ii) = x(3);
end

wr(end)
max(Ia)

figure(1);
subplot(5, 1, 1);
plot(t, Va); title('Input Signal'); ylabel('Voltage [V]'); xlabel('Time [seg]');
xlim([0, 2.5]);
subplot(5, 1, 2);
plot(t, wr); title('Angular Velocity'); ylabel('Angular Velocity [rad/seg^2]'); xlabel('Time [seg]');
xlim([0, 2.5]);
subplot(5, 1, 3);
plot(t, theta); title('Angular Position'); ylabel('Angular Posit [rad/seg]'); xlabel('Time [seg]');
xlim([0, 2.5]);
subplot(5, 1, 4);
plot(t, Ia); title('Current'); ylabel('Current [A]'); xlabel('Time [seg]');
xlim([0, 2.5]);
subplot(5, 1, 5);
plot(t, Tl); title('Torque Load'); ylabel('Torque [Nm]'); xlabel('Time [seg]');
xlim([0, 2.5]);
