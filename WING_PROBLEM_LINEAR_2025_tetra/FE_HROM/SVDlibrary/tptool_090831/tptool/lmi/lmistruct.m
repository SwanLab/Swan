function lmi = lmistruct(S, n)
%LMI

% TODO: optional arg: uncertain params

I = size(S);
sizes = I(1:end-2);
R = prod(sizes);

S = reshape(S, [R I(end-1) I(end)]);

m = I(end) - n;
p = I(end-1) - n;
A = S(:, 1:n, 1:n);
B = S(:, 1:n, n+1:n+m);
C = S(:, n+1:n+p, 1:n);
D = S(:, n+1:n+p, n+1:n+m);

X = sdpvar(n, n, 'symmetric');
M = cell(1, R);
for r = 1:R
	M{r} = sdpvar(m, n, 'full');
end

% X > 0
lmi.F = set(X > 0, 'poz def');

lmi.sizes = sizes;
lmi.n = n;
lmi.m = m;
lmi.p = p;
lmi.A = A;
lmi.B = B;
lmi.C = C;
lmi.D = D;

lmi.X = X;
lmi.M = M;
