load('pend_data', 'S', 'n');

lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.5); % |x(t)| < C * exp(-0.5 t)
lmi = lmi_input(lmi, 30, 0.3);  % max 30N force when |x| < 0.3
K = lmi_solve(lmi);

save('pend_data', '-append', 'K');
