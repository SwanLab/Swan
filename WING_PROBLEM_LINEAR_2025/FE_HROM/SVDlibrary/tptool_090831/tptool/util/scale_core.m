function S = scale_core(D, v)
%SCALE_CORE Scale each element in a tensor
%	S = SCALE_CORE(D, v)
%
%	D - multidimensional data
%	v - scale vector or matrix
%
%	S - each element in D multiplied by appropriate element in v

% TODO: comment
n = numel(v);
r = numel(D)/n;
S = reshape(D, [r n]);
for i = 1:r
	S(i,:) = S(i,:) .* v(:)';
end
S = reshape(S, size(D));
