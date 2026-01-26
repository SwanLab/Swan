function G = tp_cl(S, Asize, K)
%TP_CL C losed loop TP system
%	G = TP_CL(S, Asize, K)
%	
%	S     - core tensor of a TP model (A, B)
%	Asize - size of matrix A in the system matrix (S = [A B; C D])
%	K     - state feedback controller
%
%	G     - cl system

I = size(S);
Bsize = I(end) - Asize;
Iomega = I(1:end-2);
R = prod(Iomega);

S = reshape(S, [R I(end-1) I(end)]);
K = reshape(K, [R Bsize Asize]);

A = S(:, 1:Asize, 1:Asize);
B = S(:, 1:Asize, Asize+1:Asize+Bsize);

G = zeros(R,Asize,Asize);
for r = 1:R
	Ar = reshape(A(r,:,:), [Asize Asize]);
	Br = reshape(B(r,:,:), [Asize Bsize]);
	Kr = reshape(K(r,:,:), [Bsize Asize]);
	G(r,:,:) = Ar-Br*Kr;
end
G = reshape(G, [Iomega Asize Asize]);
