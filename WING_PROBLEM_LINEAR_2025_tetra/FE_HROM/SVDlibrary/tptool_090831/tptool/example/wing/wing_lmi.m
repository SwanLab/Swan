load('wing_data', 'S', 'n');

% State feedback TP controller design
lmi = lmistruct(S, n);
%lmi = lmi_asym(lmi);
lmi = lmi_asym_decay(lmi, 0.1);
K = lmi_solve(lmi);

save('wing_data', '-append', 'K');
