function [t, wr]= motorModel()
clear all; close all; clc;

%% State Space Model
%
%   [xp] = A*[x] + B*[u]
%   [y]  = C*[x] + D*[u]
%

ia(1)    = 0;
wr(1)    = 0; 
theta(1) = 0;

x   = [ia(1) wr(1) theta(1)]';
Xop = [0 0 0]';

Ra= 2.27; Laa= 0.0047; J= 0.00233; Bm= 0.00131;
Ki= 0.25; Km= 0.25;

A= [ -Ra/Laa  -Km/Laa  0  ;
      Ki/J    -Bm/J    0  ;
        0        1     0  ];

B= [  1/Laa       0   ;
        0       -1/J  ;
        0         0   ];

C= [0 1 0];

D= [0 0];

% size(A)
% size(B)
% size(C)
% size(D)

tsim= 5;
tint= 10e-7;

h= tsim/tint;

t= 0: tint: tsim;
Va= zeros(1, length(t));
tl= zeros(1, length(t));

for i=1: length(t)-1
    
    if t(i)>0.5
       Va(i)= 12; 
    else
       Va(i)= 0;
    end

    xp = A*(x-Xop) + B*[Va(i) tl(i)]';
    x  = x + xp*tint;
    Y  = C*(x-Xop) + D*[Va(i) tl(i)]';

    ii        = i+1;
    ia(ii)    = x(1);
    wr(ii)    = Y(1);
    theta(ii) = x(3);

end

figure(1)
%subplot(3, 1, 1); 
plot(t, wr); xlim([0, 1.5]); grid on;
title('Angular Velocity \omega_{R}');
xlabel('Time [seconds]'); ylabel('Angular velocity [rad/sec]');
%subplot(3, 1, 2);
figure(2)
plot(t, ia); grid on;
title('Current i_{a}');
xlabel('Time [seconds]'); ylabel('Current [A]');
%subplot(3, 1, 3);
figure(3)
plot(t, theta); grid on;
title('Angular Position \theta_{R}');
xlabel('Time [seconds]'); ylabel('Angular Position [rad]');

end