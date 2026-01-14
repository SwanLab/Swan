function M = hull_manip(U, Uhat)
% HULL_MANIP gives transform to achieve prescribed hull
%
% TODO...

N = size(U,2);
[i j] = find(Uhat > -1);
uhat = Uhat(Uhat > -1);
mlen = N^2;

% TODO: check if Uij in [0..1]

% min || Uvec * mvec - uhat ||

Uvec = zeros(numel(uhat), mlen);
for k = 1:numel(uhat)
	Uvec(k,(j(k)-1)*N+1:(j(k)-1)*N+N) = U(i(k),:);
end

H = Uvec'*Uvec;
f = -uhat'*Uvec;

% NN
tmp = cell(1,N);
for k = 1:N
	tmp{k} = -U;
end
Aineq = blkdiag(tmp{:});
bineq = zeros(numel(U), 1);

% SN
tmp = cell(1,N);
for k = 1:N
	tmp{k} = eye(N);
end
Aeq = horzcat(tmp{:});
beq = ones(N,1);

x0 = eye(N);
x0 = x0(:);

mvec = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],x0);
M = reshape(mvec, [N N]);
