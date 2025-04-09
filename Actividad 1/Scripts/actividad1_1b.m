%% Sistemas de control II - FCEFyN (UNC)

%% State variables representation of the system.
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u]

%% Design parameters.
R= 220;     % [Ohm].
L= 500e-3;  % [Hy]. 500mHy.
C= 2.2e-6;  % [F].
vin= 12;    % [V].

%% Operation point.
% Initial values x(0)
Il(1)  = 0;        % Current.
Vrl(1) = 0;        % Voltage through R. ??
y(1)   = 0;        % Output.

Xop    = [0, 0]';
x      = [Il(1) Vrl(1)]'; % ??

%% State Space.
A= [-R/L -1/L; 1/C 0];   % State matrix.
B= [1/L; 0];             % Input matrix.
c= [R 0];                % Output matrix.
D= 0;                    % Direct transmission matrix.

%% Obtain Transfer Function.
[num, den]= ss2tf(A, B, c, D);
G= tf(num, den)     % Build system transfer function. 
poles= pole(G)


%% Obtain natural frequency (wn), damping ratio (psita).
% Obtain these parameters to calculate integration time.
[wn, zeta]= damp(G);
damp(G)

% 0<zeta<1 underdamping system.

%% Calculate integration time (t_int), simulation time (t_sim).
% t_int calculation based wd.
% Damped natural frequency (wd) calculation.
wd= wn(1)*sqrt(1-(zeta(1))^2);
tint= ((2*pi)/(wd))/100;
tint= 0.2e-6
%10e-6

% t_sim calculation based complex pole real part. 
tsim= log(0.05)/(real(poles(1)));
tsim= 20e-3
% Define simulation step
h= tsim/tint

%% Define arrays for integration process 
% t= 0: h: tsim;
t= linspace(0, tsim, h);
u= linspace(0, 0, h);
hlf= length(t)/2; 
swTime= 0;

%% Integration
for i= 1: h-1
    swTime = swTime+1;
    if (swTime >= hlf)
        swTime = 0;
        vin = vin*(-1);
    end

    u(i) = vin;  

    xp= A*(x-Xop)+B*u(i);
    x= x+ xp*tint;
    Y= c*(x-Xop)+D*u(i);

    ii= i+1;
    y(ii)= Y(1);
    Il(ii)= x(1);
    Vrl(ii)= x(2);
end


%% Plot system response
figure(1)
plot(t, y)
ylim([-9, 5])
grid;