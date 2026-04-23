clear
close all
u3     = 0.05;
Df     =  0.2;
Dc = 0.00;
Kelas  = 1e6;
Kcoh   = 1e8;

A = -1/(Df-Dc);
B = 1 + Kelas/Kcoh;
C = -u3 * Kelas/Kcoh + Dc / (Df-Dc);

coefficients = [A, B, C];
u2 = roots(coefficients);
jump = u2;
d =(u2)/(Df); d = min(1,max(0,d));
F = jump .* (1-d) .* Kcoh;
E =0.5 * Kelas*(u3-u2).^2 + 0.5*(1-d).*(u2).^2;


% Create table with results
T = table(u2, d, F, E, ...
    'VariableNames', {'Jump', 'd', 'F', 'E'});
disp(T)


% Find index of minimum energy
[~, idx_min] = min(E);
best_u2 = u2(idx_min);
best_d  = d(idx_min);
best_F  = F(idx_min);
best_E  = E(idx_min);
fprintf('\nSolution with minimum E:\n');
fprintf('Jump = %.6f, d = %.6f, F = %.6f, E1 = %.6f\n', ...
    best_u2, best_d, best_F, best_E);

