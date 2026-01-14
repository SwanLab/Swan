function lmi = lmi_input(lmi, umax, phi)
%LMI_INPUT

% reference:
%	K. Tanaka, H. O. Wang
%	Fuzzy Control Systems Design and Analysis, (2002)
%	page 66,69 (theorem 11,13)

% TODO: known x(0), output constraint

R = size(lmi.A, 1);
X = lmi.X;
M = lmi.M;

% constraints on the control value

% phi^2 I < X
lmi.F = lmi.F + set(phi^2 * eye(lmi.n) < X, 'phi^2 I < X');

% [X, Mr'; Mr, mu^2 I] > 0
for r = 1:R
	lmi.F = lmi.F + set([X M{r}'; M{r} umax^2*eye(lmi.m)] > 0, sprintf('type3 lmi %d', r));
end
