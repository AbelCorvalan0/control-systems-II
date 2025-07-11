%% Actividad 2 - Item 3

m = 0.1;
F = 0.1;
l = 1;
g = 9.8;
M = 1.5;
%% Build arrays.

A = [ 0     1           0        0 ;   % Displacement.
      0   -F/M       -m*g/M      0 ;   % Displacement derivated.
      0     0           0        1 ;   % Angle.
      0  F/(l*M)  g*(M+m)/(l*M)  0 ];  % Angle derivated.

B = [     0     ;
         1/M    ;
          0     ; 
      -1/(l*M)  ];

C = [ 1 0 0 0 ;    % Take displacement and angle as system output.
      0 0 1 0 ];

D = [0 0]';


