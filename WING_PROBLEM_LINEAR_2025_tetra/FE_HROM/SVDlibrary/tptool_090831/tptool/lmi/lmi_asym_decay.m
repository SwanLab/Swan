function lmi = lmi_asym_decay(lmi, alpha)
%LMI_ASYM_DECAY

% reference:
%	K. Tanaka, H. O. Wang
%	Fuzzy Control Systems Design and Analysis, (2002)
%	page 62


R = size(lmi.A, 1);
A = lmi.A;
B = lmi.B;
X = lmi.X;
M = lmi.M;
n = lmi.n;
m = lmi.m;

% X*Ar' + Ar*X - Br*Mr - Mr'*Br' + 2*alpha*X < 0
for r = 1:R
	Ar = reshape(A(r,:,:), [n n]);
	Br = reshape(B(r,:,:), [n m]);
	lmi.F = lmi.F + set(X*Ar' + Ar*X - Br*M{r} - M{r}'*Br' + 2*alpha*X < 0, sprintf('type1 lmi %d', r));
end

% X*Ar' + Ar*X + X*As' + As*X - Br*Ms - Ms'*Br' - Bs*Mr - Mr'*Bs' + 4*alpha*X <= 0
for r = 1:R
	for s = r+1:R
		Ar = reshape(A(r,:,:), [n n]);
		As = reshape(A(s,:,:), [n n]);
		Br = reshape(B(r,:,:), [n m]);
		Bs = reshape(B(s,:,:), [n m]);
		lmi.F = lmi.F + set(X*Ar' + Ar*X + X*As' + As*X - Br*M{s} - M{s}'*Br' - Bs*M{r} - M{r}'*Bs' + 4*alpha*X <= 0, sprintf('type2 lmi %d', r));
	end
end
