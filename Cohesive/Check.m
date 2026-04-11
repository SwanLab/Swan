clear
clc
close all

u1     = 0.4;
Kcoh   = 1e8;
Dc     = 0.001;
Df     = 0.1;
Kelas  = 1e6;

A = Kcoh;
B = Kelas*(Df-Dc) - Kcoh*(2*u1-Df);
C = -Kcoh*(Df*u1-u1^2);

coefficients = [A, B, C];
roots = roots(coefficients);

u2 = roots(roots > 0);

jump = u1-u2;
d =(u1-u2-Dc)/(Df-Dc);
