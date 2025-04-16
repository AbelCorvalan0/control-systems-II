clc; clear all;
close all;

syms s t
syms R L C 
syms I ve vr

eq1= vr== I*R;
eq2= ve== s*L*I + I*(1/s*C) + I*R;

sol    = solve([eq1, eq2], [vr, ve]);
FT     = simplify(sol.vr/sol.ve)
[N, D] = numden(FT);

n= (R/L)*s;
d= (s^2) + s*(R/L) + (C/L);

disp('Función transferencia')
G= n/d

yt= ilaplace((1/s)*(FT));

%[N, D]= numden(FT)
%latex_str= latex(ilaplace(FT));
latex_str= latex(yt);

file= fopen('RLC equation.md', 'w');
fprintf(file, '$$\\begin{equation}\n%s\n\\end{equation}\n$$', latex_str);
fclose(file);

R= 220; L= 500e-3; C=2.2e-6;
tsim = 20e-3;
h    = 0.5e-3;
e    = 0: h: tsim;
%x= subs(yt, [R L C], [220 500e-3 2.2e-6])
%plot(t, yt)


































































