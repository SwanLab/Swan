function n = NormalsBoundary_old(Xe,ndim)
% Given the   matrix (Xe) of global coordinates of the position of nodes of boundary element "e", this function
% computes  the normal to all these elements
% Vectorized format
% ------------------
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 8-MAy-2017
if nargin == 0
    load('tmp3.mat')
end
nnodeEb = size(Xe,2) ;
nelemB = size(Xe,1)/ndim ;


n = zeros(nelemB*ndim,1) ;
% Vector t1
t1 = Xe(:,2)-Xe(:,1) ;
t1 = NormalizeVector(t1,ndim) ; % Normalization
if ndim==2
    ind2 = 2:ndim:ndim*nelemB ;
    ind1 = 1:ndim:ndim*nelemB ;
    n=[-t1(ind2);t1(ind1)];
elseif ndim == 3
    t2 = Xe(:,3)-Xe(:,1) ;    t2 =NormalizeVector(t2,ndim) ;
    %% Cross-product
    t1t2 = CrossProductVect(t1,t2,ndim) ;
    n =  NormalizeVector(t1t2,ndim)  ;
end
