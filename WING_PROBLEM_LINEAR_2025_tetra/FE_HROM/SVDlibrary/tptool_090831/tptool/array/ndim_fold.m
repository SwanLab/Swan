function T = ndim_fold(M, dim, siz)
%NDIM_FOLD Restore a matrix into a multidimensional array with the given size
%	T = NDIM_FOLD(M, dim, siz)
%	
%	M   - matrix
%	dim - the columns of M will turn into this dimension of T
%	siz - size(T) (prod(siz)==length(M), siz(dim)==size(M,1) should hold)
%
%	T   - multidimensional array containing the column vectors of M along dim dimension
%
%	eg. ndim_fold([2*ones(3) 3*ones(3)], 1, [3 3 2])
%
%	See also NDIM_UNFOLD.

%assert(prod(siz)==length(M))

ndim = length(siz);

if ndim == 2
	% T is 1D or 2D
	if dim == 2
		T = M';
	else
		T = M;
	end
else
	% TODO: better wshift
	
	% restore M into T
	new_size = wshift('1D', siz, dim-1);
	new_size(1) = size(M, 1);
	T = reshape(M, new_size);
	
	% rotate T into the right position
	dim_order = wshift('1D', [1:length(siz)], -(dim-1));
	T = permute(T, dim_order);
end
