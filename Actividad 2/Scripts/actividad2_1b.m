clear all; 
close all; clc;

%% Import data table.
% Read table.
data= readtable('../Resources/Curvas_Medidas_Motor_2025_v.xlsx');
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t), capacitor voltage(Vcap) and current(I) array.
% te   = dataT(:, 1);
% wr  = dataT(:, 2);
% I   = dataT(:, 3);
% Va  = dataT(:, 4);
% tl  = dataT(:, 5);
% disp('Lenght of tl vector')
% length(te)
% disp('Maximium value time.')
% max(te)
% 
% disp('Time simulation: ')
% tint= te(2)-te(1)

%% Define system parameters.
%
Ra  = 2.258930051299405;
La  = 0.005026901184834716;
Km  = 0.2500481104997174;
Ki  = 0.25965987053759737;
Jm  = 0.0028472626983113334;
Bm  = 0.0014165170369840668;


%% Define State Space Model.

% State Matrix.
A= [-Ra/La -Km/La  0;
     Ki/Jm -Bm/Jm  0;
       0      1    0];

% Input matrix.
B= [ 1/La ;
      0   ;
      0   ];

% Output Matrix.
C= [0 0 1];

D= [0 0];

%sys= ss(A, B, C, D)

%% Controllability verification.
% Controllability matrix calculation.
%Co= ctrb(sys)
% It allows us if the system can be carried brought into a any state
% in a limited time.

% Range of controllability matrix.
%rank_Co= rank(Co);

%if rank(Co) == size(A, 1)
%    disp('Sistema controlable');
%else
%    error('Sistema no controlable');
%end   

%% Define LQR Controller.
% Build Q matrix (weight of state error).
Q= diag([0.1, 0.1, 1, 100]);
% Build R matrix (weight of control signal).
% One input signal.
R= 1; 

% Eigen values from A.
% eigenV= eig(A)  % display eigenvalues.
% We've two negative real part eigenvalues
% stable system.
% Besides we've third eigenvalue which is equal to zero (0).
% which makes the system marginally stable.

%% Expanded matrix.
% We need to expand system to add error integral as state variable.
%
% |  A  0  |   | x(t) |   | B |        | 0 | 
% |        | * |      | + |   | u(t) + |   | r(t)
% | -C  0  |   | e(t) |   | 0 |        | 1 | 

% In this case we've size(A), size(B) and size(C) the following results:
size(A)
% size(A) = 3x3.
% We need to add: 
% 3x1 matrix of zeros.
size(B)
% size(B) = 3x2.
% 1x1 escalar because we want simulate the system only with 
% respect to u(t).
size(C)
% size(C) = 1x3.
% 1x1 escalar (zero).

Aa= [  A  zeros(3,1) ;
      -C      0      ];

Ba= [  B(:,1)  ; 
          0    ];

Ca= [ C 0 ]

%% Control signal u(t).
% The control signal is determinated by:
% u(t)= -K x(t)
% u(t)= -(R^(-1) B^(T) P) x(t)
% It's defined lambda(t).
% lambda(t)= P(t) x(t).
% P allows calculate the optimium controller gain.

%% Hamiltonian calculation.
% P is a positive simetric matrix which solves
% the Ricatti's equation.
%
%     |  A   -B R^(-1) B^(T)  |
% H = |                       |
%     | -Q        -A^(T)      |

% Calculate of 
R_inv= inv(R);

H= [  Aa -Ba*R_inv*(Ba');
      -Q       -Aa'     ];

%% Obtain eigenvectors, eigenvalues of H matrix.
% Eigen decomposition.
[V, D]= eig(H);

% Returns diagonal matrix D of eigenvalues and matrix V 
% whose columns are the corresponding right eigenvectors, 
% so that A*V = V*D.

% Get eigen values.
eigenV= diag(D);

% Analyze stable eigenvectors.
stableEv  = real(eigenV)<0;
V_estable = V(:, stableEv);

% Separate components.
size(V)
% 8x8 vector.
MX1= V_estable(1:4, :);
MX2= V_estable(5:8, :);

MX1_inv= inv(MX1);

% Obtain P matrix.
P= real(MX2*MX1_inv);

% Calculate Ka.
Ka= R_inv*(Ba')*P;
% What's that for?
% In Ka we've the following structure:
% Ka= [ K  KI ];
size_Ka= size(Ka)
%Ka
% Build K.
K= Ka(1:3);
% Build KI.
KI= Ka(size_Ka(2));

%% Verifying stability of closed-loop system.
eig(Aa-Ba*Ka)
% Stable system. It may have a slight oscillation.

% Ask for this doubt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to calculating integration time, 
% which value I must take?

%% Build vectors for simulation.
% Build vectors in function of length of "t" vector 
% which determinates the step of simulation.
h = 1e-4;
simTime = 50;
t = 0:h:(simTime-h);

%% Generate reference (theta).
% This reference signal toggle between pi/2 and -pi/2.
% Inicialización de ref (señal que cambia cada 5 segundos entre +pi/2 y -pi/2)
ref = pi/2 * ones(size(t)); % Begin with pi/2 value.
Tl  = zeros(1, length(t));

%% Initial conditions.
ia(1)    = 0;
wr(1)    = 0;
theta(1) = 0;

stateVector = [ ia(1)  wr(1)  theta(1) ]';

Xop      = [0 0 0]';

x= [ ia(1) wr(1) theta(1) ]';

zeta(1)  = 0;
integ(1) = zeta(1);

%% Observator building.
% Duallity system.
% Not augment system implemented.
% Open-loop.

Ao = A'
Bo = C'
Co = B'

Qo = diag([1, 1, 0.1])
Ro = 10;

% Ro disminuiir y agregar otra variable

Ko = lqr(Ao, Bo, Qo, Ro)

obsStateVector = [  ia(1)  wr(1)  theta(1)  ]';
xObs = [ 0 0 0 ]';

% Ca(1:3) toma estados del vector, EXCLUYE DERIVADA DEL ERROR.

for i = 1:(simTime/h)
    
    % Build torque signal (Tl)
    if t(i) >= 0.701 && t(i) <= 1.0
        Tl(i) = 0.12; % Amplitude.
    else
        Tl(i) = 0;    % Amplitude out of work time.
    end

    % Toggle reference (ref) each 5 seconds 
    % (alternate between +pi/2 y -pi/2).
    if mod(floor(t(i)), 10) < 5  % Toggle each 5 seconds.
        ref(i) =  pi/2;
    else
        ref(i) = -pi/2;
    end
    
    % Error derivated.
    zetaP   = ref(i) - Ca(1:3)*stateVector - Ca(4)*integ;
    % Integral of error.
    zeta(i) = integ + zetaP * h; 
    % Control signal.
    u(i)= -K*stateVector - KI*zeta(i);
    % Load new state values.
    ia(i)    = x(1);
    wr(i)    = x(2);
    theta(i) = x(3); 
    % Calculate of states.
    x1P = -Ra*x(1)/La - Km*x(2)/La + u(i)/La;
    x2P =  Ki*x(1)/Jm - Bm*x(2)/Jm - Tl(i)/Jm;
    x3P =  x(2);
    
    xP  = [ x1P  x2P  x3P ]'; 
    x   = x + h*xP;

    iaO(i)    = xObs(1);
    wrO(i)    = xObs(2);
    thetaO(i) = xObs(3);

    y0(i) = C*obsStateVector;
    y(i)  = Ca(1:3)*stateVector + Ca(4)*integ;
    
    xTP  = A*xObs + B*u(i) + Ko*(y(:, i) - y0(:, i));
    xObs = xObs + xTP*h;

    stateVector    = [ ia(i)  wr(i)  theta(i) ]';
    integ          = zeta(i);
    obsStateVector = [ iaO(i)  wrO(i) thetaO(i) ]';

end

%% 
% La corriente observada es:
figure(1)
plot(t, iaO,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('Corriente [A]')
title('Corriente de armadura i_a observada')
grid
%% 
% Ahora, para comparar con la real obtenida anteriormente:
figure(2)
plot(t, iaO,'LineWidth', 1.5); hold on;
plot(t, ia, 'LineWidth', 1.5); xlim([0, 15]);
xlabel('Time (seconds)'); ylabel('Current (Ampere)');
title('Armature current i_a')
legend('Observada', 'Real')
grid on
hold off
%% 
% Es posible apreciar un error en la aproximación realizada mediante el observador. 
% A continuación se muestra la evolución de ese error.

% obsError(1) = 0;
% for i = 1:(simTime/h)
%     obsError(i) = (realIa.ia(i)-iaO(i));
% end
% plot(t,obsError,'LineWidth',1.5)
% xlabel('Tiempo [seg]')
% ylabel('Error [A]')
% title('Error de observación')
% grid

% figure(1);
% subplot(4, 1, 1); plot(t, ref, '--', 'Color', [1 0.5 0]);
% hold on; plot(t, theta, 'b', 'LineWidth', 1.5); 
% xlim([0, 20]); grid on;
% title('Reference \theta(t)'); 
% xlabel('Time (Seconds)'); ylabel('Angular position (rad)')
% 
% subplot(4, 1, 2); plot(t, u, 'LineWidth', 1.5);
% xlim([0, 20]); ylim([-2.3, 2.3]); grid on;
% title('Control signal u(t)');
% xlabel('Time (Seconds)'); ylabel('Voltage (Volt)');
% 
% subplot(4, 1, 3); plot(t, Tl, 'LineWidth', 1.5);
% xlim([0, 20]); ylim([0, 0.15]); grid on;
% title('Torque T_{l}(t)');
% xlabel('Time (Seconds)'); ylabel('Torque (Nm)');
% 
% subplot(4, 1, 4); plot(t, ia, 'LineWidth', 1.5);
% xlim([0, 20]); ylim([-0.7, 0.7]); grid on;
% title('Current i_{a}(t)');
% xlabel('Time (Seconds)'); ylabel('Current (Ampere)');