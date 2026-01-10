function T = tprod(S, U)
%TPROD Tensor product of an ndim-array and multiple matrices
%	T = TPROD(S, U)
%	
%	S  - tensor (multidimensional array)
%	U  - cell containing matrices for each dimension of S
%	     if U{n}==[] or length(U)<n then no multiplication will be done
%	     in the nth dimension ie. it is the same as U{n}==eye(size(S,n))
%	
%	T  - result of the product
%
%	eg. tprod(ones(2,3,4), {ones(5,2), [], ones(6,4)})
%
%	See also TPROD1.

% TODO: longer U ->  expand

T = S;
siz = size(S);
for i = 1:length(U)
	if ~isempty(U{i})
		siz(i) = size(U{i},1);
		H = ndim_unfold(T, i);
		T = ndim_fold(U{i}*H, i, siz);
	end
end
