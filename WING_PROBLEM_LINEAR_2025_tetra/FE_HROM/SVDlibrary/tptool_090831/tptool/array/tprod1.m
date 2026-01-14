function T = tprod1(S, U, n)
%TPROD1 Tensor Product of a tensor and a matrix
%	T = TPROD1(S, U, n)
%	
%	S  - tensor (multidimensional array)
%	U  - matrix compatible with the nth size of S ie. size(S,n)==size(U,1)
%	n  - execute tprod in this dimension
%	
%	T  - result of the product
%
%	eg. tprod1(ones(2,3,4), [1 0 0; 0 1 0], 2)
%
%	See also TPROD.

% TODO: n > dimensions of S -> ndim_expand

siz = size(S);
siz(n) = size(U,1);
H = ndim_unfold(S, n);
T = ndim_fold(U*H, n, siz);
