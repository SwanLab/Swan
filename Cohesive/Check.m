clear
close all
u3 = 0.21;
Df     = 0.2;
Kelas  = 2;
Kcoh   = 4;

A = -1/Df;
B = 1 + Kelas/Kcoh;
C = -u3 * Kelas/Kcoh;
coefficients = [A, B, C];
u2 = roots(coefficients);

jump = u2;
d =(u2)/(Df);
d = min(1,max(0,d))

F = jump .* (1-d) .* Kcoh;

E =0.5 * Kelas*(u3-u2).^2 + 0.5*(1-d).*(u2).^2


