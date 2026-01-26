load('tora_data', 'S', 'n');

% State feedback TP controller design
lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
K = lmi_solve(lmi);

save('tora_data', '-append', 'K');
