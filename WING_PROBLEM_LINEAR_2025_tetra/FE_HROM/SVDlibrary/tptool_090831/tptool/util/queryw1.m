function W = queryw1(U, domain, p)
%QUERYW1 Query weight function values at given parameter point
%	W = QUERYW1(U, domain, p)
%	
%	U        - (discrete) weight function data of the TP model
%	domain   - intervals for each dimension
%	p        - parameter point where we query the weight functions
%
%	W        - linear interpollated weights

% TODO: method param to interp1 (spline,..)

if size(domain,1) ~= length(U)
	error 'size(domain,1) == length(U) must be true'
end

W = cell(1,length(U));
for i = 1:length(U)
	x = linspace(domain(i,1), domain(i,2), size(U{i},1));
	% TODO: extrapolation
	if p(i) < x(1)
		W{i} = U{i}(1,:);
	elseif p(i) > x(end)
		W{i} = U{i}(end,:);
	else
		W{i} = interp1(x, U{i}, p(i));
	end
end
