function S = shift_core(D, v)
%SHIFT_CORE Shift the elements of a tensor
%	S = SHIFT_CORE(D, v)
%
%	D - multidimensional array
%	v - shift matrix
%
%	S - shifted array

% TODO: comment
n = numel(v);
r = numel(D)/n;
S = reshape(D, [r n]);
for i = 1:r
	S(i,:) = S(i,:) + v(:)';
end
S = reshape(S, size(D));
