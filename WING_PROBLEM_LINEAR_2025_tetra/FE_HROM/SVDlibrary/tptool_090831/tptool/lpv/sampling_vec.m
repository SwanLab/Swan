function D = sampling_vec(f, domain, gridsize)
%SAMPLING_VEC Sampling a vector-vector function
%	D = SAMPLING_VEC(f, domain, gridsize)
%	
%	f        - function handle with a real vector argument with N elements
%	domain   - [min1 max1;... minN maxN] intervals for each element
%	gridsize - number of sampling grid points for each element
%	
%	D        - (N+1)-dimensional array of the sampled data
%
%	eg.  sampling_vec(@(x) [1 x(1) x(2)], [-2 2; 0 3], [3 5])
%
%	See also SAMPLING_FUNC

% TODO: now it works if f returns multidim array..

tmp = f(domain(:,1));
SIZout = size(tmp);
Nout = numel(tmp);
Nin = size(domain, 1);
if length(gridsize) == 1
	gridsize = gridsize * ones(Nin,1);
else
	if not(length(gridsize) == Nin)
		error 'length of "gridsize" and "domain" must be the same.'
	end
	gridsize = gridsize(:);
end

if Nin > 0
	% allocation of sample array
	siz = prod(gridsize);
	D = zeros(siz * Nout,1);
	a = domain(:,1);
	b = domain(:,2);
	step = (b - a) ./ (gridsize - 1);

	% sampling
	z = ones(Nin,1);
	z(1) = 0;
	for k = 1:siz
		% next gridsize point index
		z = nexti(gridsize, z);

		% arguments
		p = a + step .* (z - 1);

		% sampling the model at the gridsize point
		y = f(p);
		D((k-1) * Nout + 1 : k * Nout) = y(:);
	end

%	if Nin > 1
		% reshape to proper size (last dim is the out vec dim)
		D = shiftdim(reshape(D, [SIZout gridsize']), length(SIZout));
%	end
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
