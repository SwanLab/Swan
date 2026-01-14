function data = sampling_lpv(lpv, dep, domain, gridsize)
%SAMPLING_LPV Sampling an LPV model
%	sampled = SAMPLING_LPV(lpv, dep, domain, gridsize)
%
%	lpv      - LPV model (cell array of functions with given domain)
%	dep      - parameter dependency of each function (makes sampling more efficient)
%	domain   - sampling intervals for each parameter
%	gridsize - sampling grid size for each parameter
%
%	data     - cell array of the sampled data, not a whole ndim array
%	           (should be used together with domain and gridsize)
%
%	eg. sampling_lpv({@(x)1 @(x)x(1) @(x)x(1)+x(2)}, [0 0; 0 1; 1 1], [-2 2; 0 3], [3 5])
%
%	See also SAMPLING_VEC

% TODO: lpv docs, data docs
% TODO: data -> ndim array
% TODO: gridsize(i) == 1

%DIMENSIONS, PARAMETERS
[Sy, Sx] = size(lpv);
P = size(domain,1);
if length(gridsize) == 1
	gridsize = gridsize * ones(P,1);
else
	if not(length(gridsize) == P)
		error 'length(gridsize) == P must be true'
	end
	gridsize = gridsize(:);
end

%SAMPLING ELEMENT BY ELEMENT OF Sx * Sy MATRIX
data = cell(Sy,Sx);
if Sy == 1 || Sx == 1
	for i = 1:length(lpv)
		% index of dependent params
		depidx = squeeze(dep(i,:) > 0);
		% sampling
		data{i} = sampling_dep(lpv{i}, depidx, domain(depidx,:), gridsize(depidx));
	end
else
	for i = 1:Sy
		for j = 1:Sx
			% index of dependent params
			depidx = squeeze(dep(i,j,:) > 0);
			% sampling
			data{i,j} = sampling_dep(lpv{i,j}, depidx, domain(depidx,:), gridsize(depidx));
		end
	end
end
