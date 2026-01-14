function B = ndim_expand(A, v)
%NDIM_EXPAND expand an array in a new dimension by multiplying it with a vector
%	B = NDIM_EXPAND(A, v)
%	
%	A - multidimensional array (can be vector or matrix)
%	v - vector
%	B - add one new dimension to A as: [A.*v(1) | A.*v(2) | ... | A.*v(n)]
%
%	eg. ndim_expand(ones(2,3), [2 3 4])

s = size(A);
d = length(s);
if d==2 && s(2)==1
	v = v(:);
	if s(1)==1
		B = A*v;
	else
		B = A*v';
	end
else
	d = d + 1;
	n = length(v);
	B = v(1)*A;
	for i = 2:n
		B = cat(d, B, v(i)*A);
	end
end
