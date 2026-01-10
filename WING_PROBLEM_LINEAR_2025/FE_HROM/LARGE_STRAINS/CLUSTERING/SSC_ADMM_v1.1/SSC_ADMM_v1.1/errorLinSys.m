%--------------------------------------------------------------------------
% This function computes the maximum L2-norm error among the columns of the 
% residual of a linear system 
% Y: DxN data matrix of N data point in a D-dimensional space
% Z: NxN sparse coefficient matrix
% err: maximum L2-norm of the columns of Y-YZ 
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function err = errorLinSys(P,Z)

[R,N] = size(Z);
if (R > N) 
    E = P(:,N+1:end) * Z(N+1:end,:);
    Y = P(:,1:N);
    Y0 = Y - E;
    C = Z(1:N,:);
else
    Y = P;
    Y0 = P;
    C = Z;
end

[Yn,n] = matrixNormalize(Y0);
M = repmat(n,size(Y,1),1);
S = Yn - Y * C ./ M;
err = sqrt( max( sum( S.^2,1 ) ) );