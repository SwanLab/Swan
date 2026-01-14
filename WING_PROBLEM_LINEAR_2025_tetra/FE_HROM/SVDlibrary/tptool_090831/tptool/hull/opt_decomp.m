function [W V] = opt_decomp(U, Uhat)
% OPT_DECOMP optimal decomposition to achieve prescribed values
%
% TODO...

% N := size(U, 2) + 1 == size(Uhat, 2)
% U := [U|1] (1 is needed for SN constraint)
% unknown: M (NxN real matrix)
% W = U M
% V = pinv(M)
%
% optimization problem:
% min [(U M)ij - (Uhat)ij]^2
% with respect to
% U M 1 = 1
% (U M)ij >= 0

K = size(U,1);
N = size(U,2) + 1;

U = [U ones(K,1)];
[i j] = find(Uhat > -1);
uhat = Uhat(Uhat > -1);
mlen = N*N;

if N ~= size(Uhat,2)
	error 'size(Uhat,2) must be size(U,2)+1'
end
% TODO: check if "U x = ones" has a solution: (I - U U') 1 == 0
% TODO: check if uhatij in [0..1]

% min || Uvec * mvec - uhat ||
Uvec = zeros(numel(uhat), mlen);
for k = 1:numel(uhat)
	Uvec(k,(j(k)-1)*N+1:(j(k)-1)*N+N) = U(i(k),:);
end

H = Uvec'*Uvec;
f = -uhat'*Uvec;

% NN
Aineq = -kron(eye(N),U);
bineq = zeros(K*N,1);

% SN
%Aeq = kron(ones(1,Nh),U);
%beq = ones(K,1);
% FIX: eliminate dependent constraints (rank(Aeq|b) = N)
[u s] = svd(U','econ');
s = s(1:N,1:N);
us = u(:,1:N)*s;
Aeq = kron(ones(1,N),us');
beq = us(end,:)'; % last column of U is ones

x0 = eye(N,N);
x0 = x0(:);

mvec = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],x0);
M = reshape(mvec, [N N]);

W = U*M;
V = inv(M);    % W V == [U | 1]
V = V(:,1:N-1); % W V == U
