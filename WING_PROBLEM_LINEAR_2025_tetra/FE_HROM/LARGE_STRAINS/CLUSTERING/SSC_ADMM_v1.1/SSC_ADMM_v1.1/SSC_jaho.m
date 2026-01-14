%--------------------------------------------------------------------------
% This is the function to call the sparse optimization program, to call the 
% spectral clustering algorithm and to compute the clustering error.
% r = projection dimension, if r = 0, then no projection
% affine = use the affine constraint if true
% s = clustering ground-truth
% missrate = clustering error
% CMat = coefficient matrix obtained by SSC
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

%function [missrate,CMat] = SSC_jaho(X,r,affine,alpha,outlier,rho,s)
function [CMat,grps] = SSC_jaho(X,nclusters,DATAoffline)

if nargin == 2 
    DATAoffline = [] ; 
end
    DATAoffline = DefaultField(DATAoffline,'NumberOfIterationsSubspaceClustering',200) ; 


%if (nargin < 6)
    rho = 1;  % What is rho ?  
%end
%if (nargin < 5)
    outlier = false;
%end
%if (nargin < 4)
    alpha = 20;  % What is alpha ? 
%end
%if (nargin < 3)
    affine = false;
%end
%if (nargin < 2)
    r = 0;   % What is r ? 
%end

n = nclusters ; 
Xp = DataProjection(X,r);

if (~outlier)
     % default coefficient error threshold to stop ADMM
    % default linear system error threshold to stop ADMM
    thr = 2*10^-4; 
    maxIter = DATAoffline.NumberOfIterationsSubspaceClustering ; 
    CMat = admmLasso_mat_func(Xp,affine,alpha,thr,maxIter);
    C = CMat;
else
    CMat = admmOutlier_mat_func(Xp,affine,alpha);
    N = size(Xp,2);
    C = CMat(1:N,:);
end

CKSym = BuildAdjacency(thrC(C,rho));
grps = SpectralClustering(CKSym,n);
%missrate = Misclassification(grps,s);