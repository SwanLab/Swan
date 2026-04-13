clear
close all

u1     = 0.4;
u3     = 0;
Kcoh   = 1e8;
Dc     = 0.001;
Df     = 0.1;
Kelas  = 1e6;

% A = Kcoh;
% B = Kelas*(Df-Dc) - Kcoh*(2*u1-Df);
% C = -Kcoh*(Df*u1-u1^2);


A = -1/(Df-Dc);
B = 1+2*u1/(Df-Dc) + Kelas/Kcoh;
C = Dc/(Df-Dc) - u1*(u1+Dc+1)/(Df-Dc)+u3*Kelas/Kcoh;

coefficients = [A, B, C];
u2 = roots(coefficients);

jump = u1-u2;
d =(u1-u2-Dc)/(Df-Dc);

F = jump .* (1-d) .* Kcoh;
