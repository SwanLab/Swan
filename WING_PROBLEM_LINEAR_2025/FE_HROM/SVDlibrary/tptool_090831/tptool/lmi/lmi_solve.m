function [K P] = lmi_solve(lmi)
%LMI_SOLVE

% solve
solvesdp(lmi.F);
checkset(lmi.F)

P = inv(double(lmi.X));
R = size(lmi.A, 1);
K = zeros(R, lmi.m, lmi.n);
for r = 1:R
	K(r,:,:) = double(lmi.M{r}) * P;
end
K = reshape(K, [lmi.sizes lmi.m lmi.n]);
