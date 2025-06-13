clear all; 
close all; clc;

%% Import data table.
% Read table.
data= readtable('../Resources/Curvas_Medidas_Motor_2025_v.xlsx');
data;
% Convert table to array. Apply transpose.
dataT= table2array(data);

% Obtain time(t), capacitor voltage(Vcap) and current(I) array.
te   = dataT(:, 1);
wr  = dataT(:, 2);
I   = dataT(:, 3);
Va  = dataT(:, 4);
tl  = dataT(:, 5);
disp('Lenght of tl vector')
length(te)
disp('Maximium value time.')
max(te)

disp('Time simulation: ')
tint= te(2)-te(1)

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
B= [1/La    0 ;
     0   -1/Jm;
     0      0 ];

% Output Matrix.
C= [0 0 1];

D= [0 0];

%% Initial conditions.
Xop      = [0 0 0]';
ia(1)    = 0;
wr(1)    = 0;
theta(1) = 0;

X= [ ia(1); wr(1); theta(1) ];

sys= ss(A, B, C, D)

%% Controllability verification.
% Controllability matrix calculation.
Co= ctrb(sys)
% It allows us if the system can be carried brought into a any state
% in a limited time.

% Range of controllability matrix.
%rank_Co= rank(Co);

if rank(Co) == size(A, 1)
    disp('Sistema controlable');
else
    error('Sistema no controlable');
end   

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

Ca= [ C 0 ];

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
ts= 0.001;
tf= 10;

t  = 0: ts: tf;

%% Generate reference (theta).
% This reference signal toggle between pi/2 and -pi/2.
ref = zeros(1, length(t));

for k= 1:length(t)
    ref(k) = (-1)^(floor(t(k)/5)) * pi/2;
end
