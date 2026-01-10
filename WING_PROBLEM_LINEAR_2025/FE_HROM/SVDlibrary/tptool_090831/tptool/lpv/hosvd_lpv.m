function [S U sv tol] = hosvd_lpv(data, dep, gridsize, svtol, keep)
%HOSVD of a discretized lpv model
%	[S U]        = HOSVD_LPV(data, dep, gridsize)
%	[S U sv tol] = HOSVD_LPV(data, dep, gridsize, svtol, keep)
%
%	data     - sampled LPV model (cell array)
%	dep      - parameter dependency array
%	gridsize - sampling grid size for each parameter
%	svtol    - singular value tolerance (1e-8 by default)
%	keep     - number of kept singular values (if set then no interactive question)
%
%	S        - core tensor of the HOSVD based canonical representation
%	U        - basis functions of the canonical representation
%	sv       - singular values for each parameter
%	tol      - singular value tolerance for each parameter
%
%	See also SAMPLING_LPV


% TODO: matlab docs

%Dependencies in needed form
dep = dep2idx(dep);

%Number of parameters
P = size(dep, 3);

%Size of S matrix
[Sy, Sx] = size(data);

if nargin < 4
	svtol = 1e-8;
end
if nargin < 5 || isempty(keep)
	keep = zeros(1,P);
end

%Weighting function, singular value and tolerance allocation
U = cell(1,P);
sv = cell(1,P);
tol = cell(1,P);
gridprod = prod(gridsize);

%% Weighting functions
%TODO: i,j,k indexing..
for i = 1:P

	%Product of grid size of other dimensions
	Mprod = gridprod/gridsize(i);

	%Initial tmp and constant
	tmp = [];
	c = 0;

	%Every element of S matrix
	for j = 1:Sy
		for k = 1:Sx
			%Number of dependencies at current element
			n = sum(dep(j,k,:)>0);

			if n > 0
				% Create subtensor to reduce computation
				if dep(j,k,i) > 0
					%If dependent, layout, svd, keep column space weighting
					%functions
					lay = ndim_unfold(data{j,k}, dep(j,k,i));
					if n > 2
						%Reduce to column space
						[u, s] = svd(lay, 'econ');
						lay = u*s;
					end
					tmp = [tmp lay.*sqrt(prod(gridsize(dep(j,k,:)==0)))];
				else
					%Add constant
					c = c + sum(data{j,k}(:).^2) * prod(gridsize(dep(j,k,:)==0))/gridsize(i);
				end
			else
				% const element
				c = c + data{j,k}^2 * Mprod;
			end
		end
	end

	tmp = [tmp ones(gridsize(i),1).*sqrt(c)]; 

	% SVD based reduction of the current dimension (i)
	[Ui svi toli] = svdtrunc(tmp, svtol);
	if keep(i) <= 0
		disp('normalised singular values (and original ones):');
		for j = 1:length(svi)
			fprintf('%12.5f  (%g)\n',svi(j)/sqrt(gridprod), svi(j));
		end
		% TODO: rewrite with isstrprop
		ns = input('number of singular values to keep [all] = ');
	else
		ns = keep(i);
	end
	if ~isempty(ns) && ns < length(svi)
		Ui = Ui(:,1:ns);
		toli = svi(ns+1);
		svi = svi(1:ns);
	end
	U{i} = Ui;
	sv{i} = svi;
	tol{i} = toli;
end

%% Calculate coretensor
S = coretensor(U, data, dep);

