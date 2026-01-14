function S = querylpv(lpv, p)
%QUERYLPV Query an LPV model at given parameter point
%	S = QUERYLPV(lpv, p)
%	
%	lpv   - LPV model
%	p     - given parameter
%	
%	S     - system matrix at the given point

% TODO: convert to lti sys

S = zeros(size(lpv));
for i = 1:size(lpv,1)
	for j = 1:size(lpv,2)
		S(i,j) = lpv{i,j}(p);
	end
end
