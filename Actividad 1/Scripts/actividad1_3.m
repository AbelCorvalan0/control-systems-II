clear all; close all; clc;
%%% Obtain Parameters R, L, C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R= 127.7277;
C= 3.9434e-06;
L= 0.0027;

%% Simulate State-Space 
num= [C 0];
den= [L*C C*R 1];
s= tf('s');
G= tf(num, den)
[A, B, c, D]= tf2ss(num, den);
p= pole(G);
% Fastest dynamic pole.
p(1)
% Slower dynamic pole.
p(2)

vin= 12;
tint = log(0.95)/p(1)
tsim = log(0.05)/p(2);
T= 20e-3;
tsim= 6.5*T;
h= tsim/tint;
t= linspace(0, tsim, h);
u= linspace(0, 0, h);
var=0;

% Hacer variable las operaciones que se repiten
% por ejemplo round(T/tint)
h10m= round(T/tint);

for i=1: h-1
    if t(i)<= 0.03
        u(i)= 0;
        %rlcModel();
    elseif t(i)> T && var<=(h10m)
        u(i)= 12;
        var= var+1;
        %rlcModel();
    elseif t(i)>= T && var>(h10m)
        u(i)= -12;
        var= var+1;
        %rlcModel();
    end

    if var>(2*h10m)
        var=0;
        u(i)= 0;
        %rlcModel();
    end
    rlcModel(A, B, c, D, u);
end

% for i=1: h-1
%     if t(i)<= T/2
%         u(i)= 0;
%     elseif t(i)> T && var<=(h10m)
%         u(i)= 12;
%         var= var+1;
%     elseif t(i)> T && var=(2*h10m)
%         u(i)= -12;
%         var= var+1;
%     elseif t(i)< 4*T && var<=(3*h10m)
%         u(i)= 12;
%         var= var+1;
%     elseif t(i)<= 5*T && var<=(4*h10m)
%         u(i)= -12;
%         var= var+1;
%     elseif t(i)<= 6*T && var<=(5*h10m)
%         u(i)= 12;
%         var= 0;
%     end
% end

plot(t, u)
%%%%% Generate .csv file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%csvwrite('simulation1.csv', data);
% filename= 'Data Generated/simulation1.csv';
% 
% headers= {'Time [Seconds]', 'Amplitude [I]'};
% fi= fopen(filename, 'w');
% fprintf(fi, '%s,%s\n', headers{:});
% fclose(fi);
% writematrix(data, filename, 'WriteMode', 'append');