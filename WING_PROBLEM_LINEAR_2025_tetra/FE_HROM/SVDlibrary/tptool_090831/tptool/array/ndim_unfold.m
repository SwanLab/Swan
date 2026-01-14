function M = ndim_unfold(T, dim)
%NDIM_UNFOLD Layout an array in a given dimension to a matrix
%	M = NDIM_UNFOLD(T, dim)
%	
%	T   - multidimensional array
%	dim - dimension which becomes the column of the resulting M matrix
%
%	M   - matrix with every vector of T along the dim dimension
%
%	eg. ndim_unfold(cat(3,ones(3,4),zeros(3,4)), 2)
%
%	See also NDIM_FOLD.

ndim = ndims(T);
if ndim == 2
	% T is 2D
	if dim == 2
		M = T';
	else
		M = T;
	end
else
	% TODO: better wshift
	
	dimorder = 1:ndim;
	newdimorder = wshift('1D', dimorder, dim-1);
	T = permute(T, newdimorder);
	
	% layout the rotated T
	siz = size(T, 1);
	M = reshape(T, siz, []);
end
