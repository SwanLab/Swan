function P = topsys(S, n)
%TOPSYS Converts a polytopic TP model into psys form
%	P = TOPSYS(S, n)
%
%	S  - core tensor of a polytopic TP model
%	n  - size of the A matrix (order of the system)
%
%	P  - psys polytopic model format of robust control toolbox

I = size(S);
m = I(end) - n;
p = I(end-1) - n;

Iomega = I(1:end-2);
R = prod(Iomega);

S = reshape(S, [R I(end-1) I(end)]);

N = max(n+p+1, 7);
M = n+m+2;

P = zeros([N 1+M*R]);
P(1,1) = -Inf;
P(2,1) = 1;
P(3,1) = R;
P(4,1) = n;
P(5,1) = m;
P(6,1) = p;
P(7,1) = 10;
for r = 1:R
	Sr = S(r, :, :);
	P(1,3+(r-1)*M+n+m) = n;
	P(n+p+1,3+(r-1)*M+n+m) = -Inf;
	P(1:n+p, 3+(r-1)*M:2+(r-1)*M+n+m) = reshape(Sr, n+p, n+m);
end

