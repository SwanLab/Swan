function N = pagenorm(X, p)
%PAGENORM Page-wise matrix or vector norm.
%   N = PAGENORM(X) computes the 2-norm of each page of N-D array X:
%              N(1,1,i) = norm(X(:,:,i))
%   The output N is an N-D array with the same number of pages as X.
%
%   If X has more than 3 dimensions, then PAGENORM implicitly expands the
%   extra dimensions to calculate the norm of all combinations:
%          N(1,1,i,j,k) = norm(X(:,:,i,j,k))
%
%   N = PAGENORM(X,p) returns the p-norm of X, where p is 1, 2, or Inf.
%   * If p = 1, then N is the maximum absolute column sum of each page of X.
%   * If p = 2, then N is the maximum singular value of each page of X.
%     This is equivalent to pagenorm(X).
%   * If p = Inf, then N is the maximum absolute row sum of each page of X.
%
%   N = PAGENORM(X,"fro") returns the Frobenius norm of each page of X.
%
%   If one of the first two dimensions of X has length 1, then PAGENORM
%   computes the page-wise vector norm for the N-D paged array of vectors V.
%
%   N = PAGENORM(V) computes the 2-norm of each page of N-D array V:
%              N(1,1,i) = norm(V(:,1,i)) or N(1,1,i) = norm(V(1,:,i))
%   The output N is an N-D array with the same number of pages as V.
%
%   N = PAGENORM(V,p) returns the generalized vector p-norm of each page of
%   V, where p is any positive real value or Inf.
%
%   See also norm, vecnorm, normalize.

%   Copyright 2021 The MathWorks, Inc.

arguments
    X {mustBeFloat,mustBeNonsparse}
    p = 2 % do validation later as we might need to modify p
end

sz = size(X);
pagesAreVectors = any(sz(1:2) == 1);

p = validateP(p, pagesAreVectors);

if isempty(X) % Empty pages
    N = zeros([1, 1, sz(3:end)], 'like', real(X));
elseif all(sz(1:2)==1) % Scalar pages
    N = abs(X);
elseif pagesAreVectors
    % call vecnorm
    N = vecnorm(X, p);
elseif matlab.internal.math.partialMatch(p, 'fro')
    % Frobenius norm reshapes matrix pages to long vectors and calls
    % vecnorm with p=2
    N = vecnorm(reshape(X, [sz(1)*sz(2), 1, sz(3:end)]), 2, 1);
elseif p == 2
    % call pagesvd and use largest singular value, which is s(1,1,:)
    s = pagesvd(X);
    N = reshape(s(1,1,:), [1 1 sz(3:end)]);

    % SVD returns NaN for any non-finite. If we encounter NaN in the
    % result we need to see if the corresponding page contained NaN. If
    % not, the page contains Inf and we correct the result to Inf.
    % We use SUM over the page dimensions as efficient way of checking for
    % NaNs.
    N(isnan(N) & ~isnan(max(X, [], [1 2], 'includenan'))) = Inf;
elseif p == 1
    N = max(sum(abs(X), 1), [], 'includenan');
else
    assert(p==Inf);
    N = max(sum(abs(X), 2), [], 'includenan');
end
end

function p = validateP(p, pagesAreVectors)
% Pre-parse data and modify some silently allowed values for p
if matlab.internal.math.partialMatch(p, 'inf')
    p = Inf;
end
if pagesAreVectors && matlab.internal.math.partialMatch(p, 'fro')
    p = 2;
end

% Verify supported values for p.
% For vector pages, p must be numerical real scalars > 0 including Inf
if pagesAreVectors
    if ~(isscalar(p) && isnumeric(p) && isreal(p) && p>0)
        error(message('MATLAB:pagenorm:unknownNorm'));
    end
else
    % For matrix pages, p must be either 1, 2, Inf or 'fro'
    if ~(matlab.internal.math.partialMatch(p, 'fro') || ...
            isscalar(p) && isnumeric(p) && any(p==[1 2 Inf]))
        error(message('MATLAB:pagenorm:unknownNorm'));
    end
end
end
