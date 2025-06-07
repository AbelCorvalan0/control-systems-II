clear all; close all;
clc;

color_= 'r-';

X        = -[0; 0; 0; 0];
t_etapa  = 60e-7;
thetaRef = 1;
tF       = 4;

Kp = 0.1;
Ki = 0.01;
Kd = 5;

Ts= t_etapa;

A1= ( 2*Kp*Ts + Ki*Ts^(2) + 2*Kd)/(2*Ts);
B1= (-2*Kp*Ts + Ki*Ts^(2) - 4*Kd)/(2*Ts);
C1= (Kd/Ts);

e  = zeros(floor(tF/t_etapa), 1); u= 0;
Tl = 0;

ii = 0;

for t= 0: t_etapa: tF
    ii = ii+1;
    k  = ii+2;
    if (2.2<=t) && (t<0.9)
        Tl= 0.12;
    end

    %% Control signal bounded.
    if u>2         
        u=2;
    end
    if u<-2
        u=-2;
    end
    
    X= modMotor(t_etapa, X, u, Tl);

    e(k)    = thetaRef-X(4);
    u       = u + A1*e(k) + B1*(k-1) + C1*e(k-2);
    acc(ii) = u;
    ia      = X(3);
    theta   = X(4);
end

t=0:t_etapa:tF;

figure(1)
subplot(4, 1, 1);hold on;
plot(t,theta,color_);title('\Theta_t [rad/sec]');hold on;

subplot(4, 1, 2); hold on;
plot(t, acc,color_);title('Accion de control');hold on;

subplot(4, 1, 3); hold on;
plot(t, ia, color_);title('Ia [A]');hold on;

subplot(4, 1, 4); hold on;
plot(t, torque, color_);title('T_L');hold on;