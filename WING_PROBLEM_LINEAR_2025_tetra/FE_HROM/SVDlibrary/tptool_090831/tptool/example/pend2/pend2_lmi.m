load('pend2_data', 'S', 'n');

lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.6);
lmi = lmi_input(lmi, 50, 0.1);
K = lmi_solve(lmi);

save('pend2_data', '-append', 'K');

