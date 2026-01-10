function S = coretensor(U, D, dep)
%Calculate HOSVD core tensor from orthonormal weighting functions and sampled data
%	S = CORETENSOR(U, D, dep)
%
%	U   - sampled basis functions for each parameter
%	D   - sampled LPV model
%	dep - parameter dependency of sampled LPV model
%
%	S   - core tensor so that D = S tprod_n U{n}

% TODO: matlab docs..

%% Initialization

%Dependencies in needed form
dep = dep2idx(dep);

%Number of parameters
P = size(dep, 3);

%Size of S matrix
siz = size(D);
Sy = siz(end-1);
Sx = siz(end);

%Pseudoinverse of weighting functions allocation
UP = cell(1,P);

%Coefficients allocation
coeff = cell(1, P);

%% Core tensor

%Pseudoinverse of weighting functions
for i=1:P
	UP{i} = pinv(U{i});
end

% Pack subtensors with the calculated weighting functions
for i = 1:P
	% setting const coeff
	coeff{i} = sum(UP{i}, 2);

	% looking for the related elements
	[depY, depX] = find(dep(:,:,i) > 0);

	% layout the corresponding tensor in the proper dimension
	for j = 1:length(depX)
		lay = ndim_unfold(D{depY(j),depX(j)}, dep(depY(j),depX(j),i));
		tmp = UP{i} * lay;
		D{depY(j),depX(j)} = ndim_fold(tmp, dep(depY(j),depX(j),i), size(D{depY(j),depX(j)}));
	end
end

% determine the final size of the core tensor S
size_W = zeros(1,P);
for i = 1:P
	size_W(i) = size(U{i}, 2);
end
size_S = [size_W Sy Sx];

Sconst = 1;
for i = 1:P
	Sconst = ndim_expand(Sconst, coeff{i});
end

% calculate the core tensor
S = zeros([Sy Sx size_W]);

% Every element of S
for i =1:Sy 
	for j = 1:Sx

		%If dependency
		if any(dep(i,j,:))
			Sp = D{i,j};
			d = ndims(Sp);
			if d==2 && length(Sp)==numel(Sp)
				d = 1;
				Sp = reshape(Sp, [length(Sp),1]);
			end
			perm = 1:P;
			z = d + 1;
			nz = 1;

			%Every dimension
			for p = 1:P
				if dep(i,j,p) == 0
					% subtensor D{i,j} does not depend on the parameter,
					% so we have to extend it along this dimension
					perm(p) = z;
					z = z + 1;
					Sp = ndim_expand(Sp, coeff{p});
				else
					perm(p) = nz;
					nz = nz + 1;
				end
			end
			%TODO: P == 1 case
			if P>1
				Sp = permute(Sp, perm);
			end
		else
			% S(i,j) is constant
			Sp = Sconst .* D{i,j};
		end

		S(i,j,:) = Sp(:);
	end
end

%Repermute tensor
length_size_S = length(size_S);
S = permute(S, [3:length_size_S 1 2]);
