function D = sampling_func(f, domain, gridsize)
%SAMPLING_FUNC Sampling a (multivariate) function
%	D = SAMPLING_FUNC(f, domain, grid)
%	
%	f      - function handle with a real vector argument with N elements
%	domain - [min1 max1;... minN maxN] intervals for each element
%	grid   - number of sampling grid points for each element
%	
%	D      - N-dimensional array of the sampled data
%
%   eg.:   sampling_func(@(x) x(1)+x(2), [-1 1; 0 3], [7 5])

N = size(domain, 1);

if length(gridsize) == 1
	gridsize = gridsize * ones(N,1);
else
	if not(length(gridsize) == N)
		error 'length of "gridsize" and "domain" must be the same.'
	end
	gridsize = gridsize(:);
end

if N > 0
	% allocation of sample array
	siz = prod(gridsize);
	D = zeros(siz,1);
	a = domain(:,1);
	b = domain(:,2);
	step = (b - a) ./ (gridsize - 1);

	% sampling
	z = ones(N,1);
	z(1) = 0;
	for k = 1:siz
		% next gridsize point index
		z = nexti(gridsize, z);

		% arguments
		p = a + step .* (z - 1);

		% sampling the model at the gridsize point
		D(k) = f(p);
	end

	if N > 1
		% reshape to proper size
		D = reshape(D, gridsize');
	end
else
	D = f([]);
end

% next gridsize point index
function n = nexti(gridsize, n)
i = 1;
n(i) = n(i) + 1;
while n(i) > gridsize(i)
	n(i) = 1;
	i = i + 1;
	n(i) = n(i) + 1;
end
