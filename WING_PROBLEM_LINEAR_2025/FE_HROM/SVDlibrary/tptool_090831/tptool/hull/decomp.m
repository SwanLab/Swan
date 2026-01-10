function [W V] = decomp(U, wtype)
%DECOMP Decomposing a matrix with orthonormal columns
%	[W V] = DECOMP(U, wtype)
%	
%	U     - matrix with orthogonal column vectors ( U' * U = I )
%	wtype - decomposition type (weight type) (default: 'snnn')
%	
%	W * V = U
%	W     - "weight" matrix
%	V     - "vertex" matrix
%	
%	wtype parameter controls the type of the column vectors in the
%	resulting W matrix.
%	If W is a "hull" type matrix it means that its row sums are 1 and each
%	element is positive ie.: W*ones == ones && all(all(W >= 0)), in this
%	case the polytope given by the rows of V contains the rows of U (the
%	convex hull of rows(V) contains the convex hull of rows(U)).
%
%	possible values of wtype:
%		'eye': hull type, W=eye(size(U,1)), V=U
%		'ortho': W orthogonal, W=U, V=eye(size(U,2))
%		'snnn' or 'hull': hull type (uses the least possible vertices)
%		'cno': 'snnn' and each column of W comes close to 1 (tight hull)
%		'irno': 'snnn' and each column of W comes close to 0
%		'box': perform box_decomp
%	
%	See also HOSVD, SNNN_DECOMP, BOX_DECOMP.

% TODO: fix cno,..
% TODO: different types..

if nargin <= 1
	wtype = 'close';
end


if ~strcmp(wtype,'eye') && ~strcmp(wtype,'ortho')
	% using an affine subspace for convex decompositions
	[n1 n2] = size(U);
	Ushift = mean(U);
	[Uu Su Vu] = svd(U - ones(n1,1)*Ushift, 'econ');
	sv = diag(Su);
	ns = sum(sv > 1e-5); % TODO
	Su = Su(1:ns,1:ns);
	Uu = Uu(:,1:ns);
	Vu = Vu(:,1:ns);

	if ns < n2
		U = Uu;
	end
end


% TODO: simplify
if nargout < 2
	switch wtype
		case 'close'
			W = close_decomp(U);
		case {'snnn', 'hull'}
			W = snnn_decomp(U);
		case 'cnoy'
			W = snnn_decomp(U);
			W = no(W);
		case 'inov'
			W = snnn_decomp(U);
			W = ino(W);
		case 'irno'
			W = snnn_decomp(U);
			W = rnoino(W);
		case 'cno'
			W = snnn_decomp(U);
			% W = cno(W, 0.5, 4, 20, 3, 3);
%			W = cno(W, 0.5, 5, 50, 10, 10);
			W = cno(W, 1, 5, 50, 15, 10);
		case 'eye'
			W = eye(size(U,1));
		case 'ortho'
			W = U;
		case 'box'
			W = box_decomp(U);
		otherwise
			error('unknown wtype');
	end
else
	switch wtype
		case 'close'
			[W V] = close_decomp(U);
		case {'snnn', 'hull'}
			[W V] = snnn_decomp(U);
		case 'cnoy'
			[W V1] = snnn_decomp(U);
			[W V2] = no(W);
			V = V2*V1;
		case 'inov'
			[W V1] = snnn_decomp(U);
			[W V2] = ino(W);
			V = V2*V1;
		case 'irno'
			[W V1] = snnn_decomp(U);
			[W V2] = rnoino(W);
			V = V2*V1;
		case 'cno'
			[W V1] = snnn_decomp(U);
			% [W V2] = cno(W, 0.5, 4, 20, 3, 3);
			[W V2] = cno(W, 1, 5, 50, 15, 10);
			V = V2*V1;
		case 'eye'
			W = eye(size(U,1));
			V = U;
		case 'ortho'
			W = U;
			V = eye(size(U,2));
		case 'box'
			[W V] = box_decomp(U);
		otherwise
			error('unknown wtype');
	end

	% shift back (for convex decompositions)
	if ~strcmp(wtype,'eye') && ~strcmp(wtype,'ortho') && ns < n2
		V = V*Su*Vu' + ones(size(V,1),1)*Ushift;
	end
end
