function C = pagetensorprod(A,B,dimA,dimB,ndimsA,ndimsB)
% paged tensor product function
%
% Syntax:
%   C = pagetensorprod(A,B);
%   ...
%   C = pagetensorprod(A,B,dimA,dimB,ndimsA,ndimsB);
%
% Inputs:
%
%   A, B: numeric arrays. (pagetensorprod is implemented to also work with
%   class objects A and/or B that support array operations). For positive
%   j, size(A,ndimsA+j) and size(B,ndimsB+j) must be matched unless either
%   size is 1.
%
%   dimA (optional, default = []): either a vector of distinct positive
%   integers, or [], or 'all', which translates to 1:ndims(A).
%
%   dimB (optional, default = dimA): either a vector of distinct positive
%   integers, or [], or 'all', which translates to 1:ndims(B). length(dimB)
%   must match length(dimA), and size(A,dimA) must match size(B,dimB).
%
%   ndimsA (optional, default = max([ndims(A);dimA(:)])): positive integer,
%   number of dimensions in A to include in the tensor product. ndimsA must
%   not be less than any element of dimA. (ndimsA corresponds to the
%   tensorprod functions' NumDimensionsA parameter, but unlike tensorprod,
%   ndimsA may be less than ndims(A) in pagetensorprod.)
%
%   ndimsB (optional, default = max([ndims(B);dimB(:)])): positive integer,
%   number of dimensions in B to include in the tensor product. ndimsB must
%   not be less than any element of dimB.
%
% Output:
%
%   C: numeric array (or a class object if either A or B is a nonempty
%   class object). If ndimsA >= ndims(A) and ndimsB >= ndims(B), then
%   pagetensorprod(A,B,...) is equivalent to tensorprod(A,B,...)
%   (generalized for class objects A and/or B). Otherwise, C is constructed
%   as follows:
%   * Take the inner product of dimensions dimA in A and corresponding
%   dimensions dimB in B.
%   * Take the outer product of dimensions setdiff(1:ndimsA,dimA) in A and
%   dimensions setdiff(1:ndimsB,dimB) in B.
%   * Apply the above operations elementwise in dimensions ndimsA+[1,2,...]
%   in A and corresponding dimensions ndimsB+[1,2,...] in B (with implicit
%   singleton expansion).
narginchk(2,6)
if nargin<3 || isequal(dimA,[])
    dimA = zeros(1,0);
elseif isequal(dimA,'all')
    dimA = 1:ndims(A);
else
    if ~(isnumeric(dimA) && isreal(dimA) && isvector(dimA) && ...
            all(isfinite(dimA) & dimA==fix(dimA) & dimA>0) && ...
            all(accumarray(dimA(:),1)<=1))
        error('pagetensorprod:validation', ...
            ['In pagetensorprod(A,B,dimA,...): dimA must be either a ' ...
            'vector of distinct positive integers, or [], or ''all''.'])
    end
    if ~isrow(dimA)
        dimA = dimA.';
    end
end
if nargin<4
    dimB = dimA;
elseif isequal(dimB,[])
    dimB = zeros(1,0);
elseif isequal(dimB,'all')
    dimB = 1:ndims(B);
else
    if ~(isnumeric(dimB) && isreal(dimB) && isvector(dimB) && ...
            all(isfinite(dimB) & dimB==fix(dimB) & dimB>0) && ...
            all(accumarray(dimB(:),1)<=1))
        error('pagetensorprod:validation', ...
            ['In pagetensorprod(A,B,dimA,dimB,...): dimB must be ' ...
            'either a vector of distinct positive integers, or [], ' ...
            'or ''all''.'])
    end
    if ~isrow(dimB)
        dimB = dimB.';
    end
end
if ~isequal(size(A,dimA),size(B,dimB))
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,...): dimA and dimB must be ' ...
        'length-matched and size(A,dimA) must match size(B,dimB).'])
end
if nargin<5
    ndimsA = max([ndims(A),dimA]);
elseif ~(isnumeric(ndimsA) && isreal(ndimsA) && isscalar(ndimsA) && ...
        isfinite(ndimsA) && ndimsA==fix(ndimsA) && ndimsA>0)
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,ndimsA,...): ' ...
        'ndimsA must be a positive integer.'])
end
if any(ndimsA<dimA)
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,ndimsA,...): ' ...
        'ndimsA must not be less than any element of dimA.'])
end
if nargin<6
    ndimsB = max([ndims(B),dimB]);
elseif ~(isnumeric(ndimsB) && isreal(ndimsB) && isscalar(ndimsB) && ...
        isfinite(ndimsB) && ndimsB==fix(ndimsB) && ndimsB>0)
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,ndimsA,ndimsB): ' ...
        'ndimsB must be a positive integer.'])
end
if any(ndimsB<dimB)
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,ndimsA,ndimsB): ' ...
        'ndimsB must not be less than any element of dimB.'])
end
% Get elementwise-multiplication dimensions sizeA_ in A, sizeB_ in B, and
% sizeC_ in C. (Singleton expansion is applied in sizeC_.)
sizeA_ = size(A,ndimsA+1:ndims(A));
sizeB_ = size(B,ndimsB+1:ndims(B));
sizeA_(end+1:length(sizeB_)) = 1;
sizeB_(end+1:length(sizeA_)) = 1;
tf = sizeA_==sizeB_ | sizeA_==1 | sizeB_==1;
if ~all(tf)
    error('pagetensorprod:validation', ...
        ['In pagetensorprod(A,B,dimA,dimB,ndimsA,ndimsB): ' ...
        'size(A,ndimsA+j) and size(B,ndimsB+j) must be matched unless ' ...
        'either size is 1. (j = ' mat2str(find(~tf)) ')'])
end
sizeC_ = sizeA_;
tf = sizeC_==1;
sizeC_(tf) = sizeB_(tf);
% Reorder A dimensions 1:ndimsA with uncontracted (outer-product)
% dimensions first and contrcted (inner-product) dimensions (dimA) last;
% then flatten the uncontracted and contracted dimensions.
p = 1:ndimsA;
p(dimA) = [];
sizeC = size(A,p);
p = [p,dimA,ndimsA+1:ndims(A)]; % dimension permutation
size_ = size(A);
size_(end+1:ndimsA) = 1;
size_ = size_(p); % dimension-permuted array size
A = permute(A,p); % size(A,j) is now size_(j).
size_ = [prod(size_(1:ndimsA-length(dimA))), ... % uncontracted dims
    prod(size_(ndimsA-length(dimA)+1:ndimsA)), ... % contracted dims
    size_(ndimsA+1:end)];
A = reshape(A,size_); % 3-D
% Reorder B dimensions 1:ndimsB with contracted (inner-product) dimensions
% (dimB) first and uncontrcted (outer-product) dimensions last; then
% flatten the contracted and uncontracted dimensions.
p = 1:ndimsB;
p(dimB) = [];
sizeC = [sizeC,size(B,p),sizeC_];
p = [dimB,p,ndimsB+1:ndims(B)]; % dimension permutation
size_ = size(B);
size_(end+1:ndimsB) = 1;
size_ = size_(p); % dimension-permuted array size
B = permute(B,p); % size(B,j) is now size_(j).
size_ = [prod(size_(1:length(dimB))), ... % contracted dims
    prod(size_(length(dimB)+1:ndimsB)), ... % uncontracted dims
    size_(ndimsB+1:end)];
B = reshape(B,size_); % 3-D
% Use pagemtimes to do the inner/outer product operations; expand flattened
% outer-product dimensions.
C = pagemtimes(A,B);
C = reshape(C,[sizeC,1,1]);
