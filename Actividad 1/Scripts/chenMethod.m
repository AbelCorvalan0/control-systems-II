function [X]= chenMethod(y1, y2, y3, g, stepAmplitude, delayTime, t_t1, zr)
% Chen's method (Second-Order Systems with differents poles).
% clear all; close all; clc;

K= g/stepAmplitude

y_1= (y1/stepAmplitude)
y_2= (y2/stepAmplitude)
y_3= (y3/stepAmplitude)

% Calculating k1, k2, k3.
% k1= ((y1)/K)-1;
% k2= ((y2)/K)-1;
% k3= ((y3)/K)-1;
k1= ((y_1)/K)-1;
k2= ((y_2)/K)-1;
k3= ((y_3)/K)-1;


% Calculating b, alfa1, alfa2.
b= 4*(k1^(3))*k3- 3*(k1^2)*(k2^2)- 4*(k2^3)+ (k3^2)+ 6*k1*k2*k3;
alfa1= (k1*k2+ k3- sqrt(b))/(2*((k1^2)+ k2));
alfa2= (k1*k2+ k3+ sqrt(b))/(2*((k1^2)+ k2));

% Calculating Beta.
% beta= (2*(k1^3)+ 3*k1*k2+ k3- sqrt(b))/(sqrt(b));
% Alternative
beta= (k1+alfa2)/(alfa1-alfa2);

% Calculating the estimates of the time constants T1, T2 and T3.
T1_e= -(t_t1-delayTime)/(log(alfa1))
T2_e= -(t_t1-delayTime)/(log(alfa2)) 
T3_e= beta*(T1_e- T2_e)+ T1_e
% T1_e= -(t_t1)/(log(alfa1));
% T2_e= -(t_t1)/(log(alfa2));
% T3_e= beta*(T1_e- T2_e)+ T1_e;

% Add estimates of the time constants T1, T2, T3.
ii= 1;
T1(ii)= T1_e;
T2(ii)= T2_e;
T3(ii)= T3_e;

% Enhancing estimation accuracy.
T1_e= sum(T1/length(T1));
T2_e= sum(T2/length(T2));
T3_e= sum(T3/length(T3));

% Build Transfer Function
s= tf('s');
switch zr
    case 0
        sys= (K*exp(-s*delayTime))/((T1_e*s+1)*(T2_e*s+1));
    case 1
        sys= (K*(T3_e*s+1))/((T1_e*s+1)*(T2_e*s+1));
    case 2
        sys= (K)/((T1_e*s+1)*(T2_e*s+1));
    otherwise
        disp('ERROR');
end

X= [sys];
end