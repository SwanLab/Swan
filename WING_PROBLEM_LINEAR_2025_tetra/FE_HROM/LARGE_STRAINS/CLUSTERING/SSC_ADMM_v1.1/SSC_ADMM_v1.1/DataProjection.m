%--------------------------------------------------------------------------
% This function takes the D x N data matrix with columns indicating
% different data points and project the D dimensional data into a r
% dimensional subspace using PCA.
% X: D x N matrix of N data points
% r: dimension of the PCA projection, if r = 0, then no projection
% Xp: r x N matrix of N projectred data points
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------


function Xp = DataProjection(X,r)

if (nargin < 2)
    r = 0;
end

if (r == 0)
    Xp = X;
else
    [U,~,~] = svd(X,0);
    Xp = U(:,1:r)' * X;
end
