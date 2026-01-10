function [W V] = close_decomp(U, n)
% CLOSE_DECOMP close hull decomposition
%
% TODO...

nU = size(U,2);
n = nU + 1;

mU = size(U,1);
idx = floor(linspace(1,mU,n));
Uhat = -ones(mU,n);
for i = 1:n
	Uhat(idx(i),i) = 1;
end

if n <= nU
	error 'size(W,2) cannot be less than the rank of U'
elseif n == nU + 1
	% usual case: nU + 1 vertex
	[W V] = opt_decomp(U,Uhat);
else
	% zeros in Uhat -> tangent
	%Uhat = zeros(mU,n);
	%for i = 1:n-1
	%	Uhat(idx(i):idx(i+1),i:i+1) = -ones(idx(i+1)-idx(i)+1, 2);
	%end
	%for i = 1:n
	%	Uhat(idx(i),i) = 1;
	%end
	
	% more vertex: many DoF --> heuristics
	% we add n-nU-1 columns to U so opt_decomp will factorize [U|A|1]
	% ideal W: Bezier curve weights: polynomials of (t + (1-t))^n == 1
	t = linspace(0,1,mU)';
	t1 = ones(mU,1)-t;
	Wb = zeros(mU,n);
	N = n-1;
	for i = 0:N
		Wb(:,i+1) = nchoosek(N,i)*(t.^i).*(t1.^(N-i));
	end
	% [U|c] can represent ones:
	c = ones(mU,1)-U*(U'*ones(mU,1));
	c = c/norm(c); % unit length
	Uc = [U c];
	% missing vectors to represent the ideal W:
	[Umissing s] = svd(Wb - Uc*(Uc'*Wb),'econ');
	U = [U Umissing(:,1:N-nU)];
	[W V] = opt_decomp(U,Uhat); % W V = [U | Umissing]
	V = V(:,1:nU);              % W V = U
	
	% shift faces -> tight (heuristics?)
	% .. 2d only for now
	if size(V,2) == 2
		V = [V(end-1,:); V(end,:); V; V(1,:)];
		ev = zeros(size(V))';
		for i = 1:size(V,1)-1
			v = (V(i,:)-V(i+1,:))';
			ev(:,i) = v./sqrt(v'*v);
		end
		for i = 2:size(V,1)-1
			n = V(i,:)' - ev(:,i)*(ev(:,i)'*V(i,:)');
			en = n/norm(n);
			d = U(:,1:2)*en;
			mod = n-max(d)*en;
			V(i,:) = V(i,:) - ev(:,i-1)'*(ev(:,i-1)'*mod);
			V(i+1,:) = V(i+1,:) - ev(:,i+1)'*(ev(:,i+1)'*mod);
		end
		V = V(3:end-1,:);
		% TODO:
		%W = smooth_w(U(:,1:2), V);
	end
end

