function D = sampling_dep(f, depidx, domain, gridsize)
%SAMPLING_DEP Sampling a (multivariate) function
%	D = SAMPLING_DEP(f, depidx, domain, gridsize)
%
%	f        - function of a real^P vector (f depends on a subset of its input)
%	depidx   - index of the dependent elements
%	domain   - [min1 max1;... minP maxP] intervals for each element
%	gridsize - number of sampling grid points for each dependent element
%
%	D        - P-dimensional array of the sampled data
%
%	This function is used internally by the tptoool toolbox
%
%	eg.:   sampling_dep(@(x) x(2)+x(3), [2 3] [-1 1; 0 3], [7 5])
%
%	See also SAMPLING_LPV

% internal function (used by sampling_lpv)

%DIMENSIONS
P = length(gridsize);
p = zeros(max(depidx), 1);

%SAMPLING
if P > 0
	% allocation of sample array
	siz = prod(gridsize);
	D = zeros(siz,1);
	a = domain(:,1);
	b = domain(:,2);
	step = (b - a) ./ (gridsize - 1);
	
	% sampling
	z = ones(P,1);
	z(1) = 0;
	for k = 1:siz
		% next grid point index
		z = nexti(gridsize, z);
		% argument
		p(depidx) = a + step .* (z - 1);
		% sampling
		D(k) = f(p);
	end

	if P > 1
		% reshape to proper size
		D = reshape(D, gridsize');
	end
else
	D = f(p);
end

%SUBFUNCTION: Next grid point index
function n = nexti(gridsize, n)
i = 1;
n(i) = n(i) + 1;
while n(i) > gridsize(i)
	n(i) = 1;
	i = i + 1;
	n(i) = n(i) + 1;
end
