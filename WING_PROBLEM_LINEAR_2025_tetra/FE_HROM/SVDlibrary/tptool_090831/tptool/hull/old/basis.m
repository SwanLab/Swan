function [U, V] = basis(A, wtype, ia)
%BASIS Basis of the column space of a matrix
%	[U,V] = BASIS(A)
%	[U,V] = BASIS(A, wtype)
%	[U,V] = BASIS(A, wtype, ia)
%%	[U,V] = BASIS(A, wtype, ia, tol)
%	
%	A     - matrix
%	wtype - weight function type (default: 'canonic')
%	ia    - interactive (default: 0)
%%	tol   - rank tolerance (default: 1e-6)
%	
%	U     - basis (weight functions)
%	V     - linear combination coefficients to get A from U (A = U*V)
%	
%	wtype parameter controls the weight function types, since (in TP model
%	transformation) one usually does not want orthonormal weight functions
%	(the result of a normal svd).
%	possible values of wtype:
%		'none': no transformation (weight functions: eye(n))
%		'canonic': svd orthonormal weigth functions (canonical form)
%		'snnn': non-negative sum-normalized (convex combination)
%		'cno': 'snnn' and each weight function comes close to 1 (tight hull)
%		'irno': 'snnn' and each weight function comes close to 0
%	
%	See also HOSVD, HOSVD_FULL.

% TODO: use tol
if nargin <= 1
	wtype = 'canonic';
end
if nargin <= 2
	ia = 0;
end
%if nargin <= 3
%	tol = 1e-6;
%end
switch wtype
	case 'canonic'
		[U, V] = svdbasis(A, ia);
	case 'snnn'
		[U, V] = snnnbasis(A, ia);
	case 'cnoy'
		[U, V] = nobasis(A, ia);
	case 'inov'
		[U] = snnnbasis(A, ia);
		[U, V] = ino(U);
	case 'irno'
		[U] = snnnbasis(A, ia);
		[U, V] = rnoino(U);
	case 'cno'
		[U] = snnnbasis(A, ia);
		% [U, V] = cno(U, 0.5, 4, 20, 3, 3);
		[U, V] = cno(U, 0.5, 5, 50, 10, 10);
	case 'none'
		U = eye(size(A,1));
		V = A;
	case 'load'
		% TODO: ugly load
		U = load('load_u.mat');
		V = U\A;
	otherwise
		error('unknown wtype');
end
