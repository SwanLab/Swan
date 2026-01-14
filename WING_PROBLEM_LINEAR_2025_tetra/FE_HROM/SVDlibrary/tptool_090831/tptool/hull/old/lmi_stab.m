function lmi_stab(A, alpha)
%LMI_STAB - check lyapunov stab. for polytopic sys. (alpha = decay rate)

I = size(A);
if (I(end) ~= I(end-1))
	error('system matrix should be square.');
end
Asize = I(end);
Iomega = I(1:end-2);
R = prod(Iomega);
A = reshape(A, [R I(end-1) I(end)]);

P = sdpvar(Asize, Asize, 'symmetric');
F = set(P > 0);
for r = 1:R
	Ar = reshape(A(r,:,:), [Asize Asize]);
	F = F + set(Ar'*P + P*Ar  < -alpha*eye(Asize));
end
solvesdp(F);
checkset(F)

