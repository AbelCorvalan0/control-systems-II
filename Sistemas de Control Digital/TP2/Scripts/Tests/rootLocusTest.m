clear all; clc; close all;

pkg load control

s   = tf('s');

num = [1, 6]
T   = conv([1, 2], [1, 4]);
den = conv(T, [1, 9]);

G= tf(num, den)

%% 1. Open poles and zeros.
pole(G)
zero(G)
%% 2. Real axis line segments.
%% 3. Infinite poles and zeros.
% We've nZ= 1 and nP= 3 in open-loop system.
% unbalanced number between poles and zeros.
% We've to introduce two infinite zeros.
% Every pole and zeros is needed to be connected.
% 2 infinite zeros, need to join to poles -2 and -4.

%% 4. Asymptotes
% Value.
% \sigma_{a}= ( sum(Fp) - sum(Fz) ) / (nP - nZ)
% \sigma_{a}= (-2+(-4)+(-9)-(-6))/(3-1)= -4.5
% Angle.
% \theta_{a}= (2*k+1)*180 /(nP -nZ)
% k=0 \theta_{a}= 180/(3-1) = 90°
% k=1 \theta_{a}= (2*1+1)*180/(3-1)= 270°

%% 5. Breakaway/ break in.
% sumM 1/(sigma-zi)= sumN 1/(sigma-pi).
% get solutions.
% sig1= -3.09 (take this value), sig2= -6.71-1.9i, sig3= -6.71+1.9i

rlocus(G)




