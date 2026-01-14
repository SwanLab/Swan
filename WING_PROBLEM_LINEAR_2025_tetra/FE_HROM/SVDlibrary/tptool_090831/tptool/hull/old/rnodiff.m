function D = rnodiff(U, c)

U0 = U(:,1:end-1);
[n,r] = size(U0);
s = sum(U0')';
U1 = U0 + (c-1)/r*s*ones(1,r);
U2 = U1 - ones(n,1)*min(U1);
U3 = U2./(ones(n,1)*max(U2));
U4 = U3./max(sum(U3'));
D = 1-min(sum(U4'))-max(U4(:,1));
