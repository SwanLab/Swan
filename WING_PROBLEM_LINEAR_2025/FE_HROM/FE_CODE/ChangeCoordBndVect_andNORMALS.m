function [Se normals]= ChangeCoordBndVect_andNORMALS(Xe,ndim)
% Given the   matrix (Xe) of global coordinates of the position of nodes of boundary element "e", this function
% expresses these positions in a coordinate system intrinsic to the
% element (Se).
% Vectorized format
% ------------------
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 27-Oct-2015
if nargin == 0
    load('tmp0.mat')
 %   ndim = 3; 
end
nnodeEb = size(Xe,2) ;
ndimB = ndim-1;
nelemB = size(Xe,1)/ndim ;  


% Rotation matrix 
R = zeros(nelemB*ndim,ndim) ; 
normals =  zeros(nelemB*ndim,1) ; 
% Vector t1
t1 = Xe(:,2)-Xe(:,1) ;   
t1 = NormalizeVector(t1,ndim) ; % Normalization
if ndimB==1
    ind2 = 2:ndim:ndim*nelemB ; 
    ind1 = 1:ndim:ndim*nelemB ; 
    R(:,1) = t1 ;  % R([ind1 ;ind2],2)=[-t1(ind2);t1(ind1)];
    normals(ind1) = t1(ind2) ; 
     normals(ind2) = [-t1(ind1)] ; 
     R(:,2) = normals ; 
elseif ndimB == 2
    R(:,1) = t1 ;
    t2 = Xe(:,3)-Xe(:,1) ;    t2 =NormalizeVector(t2,ndim) ;
    %% Cross-product 
     t1t2 = CrossProductVect(t1,t2,ndim) ;     
    R(:,3) =  NormalizeVector(t1t2,ndim)  ; 
    R(:,2) = CrossProductVect(R(:,3),R(:,1),ndim) ;    
    normals = R(:,3) ; 
end
% Transpose of R 
%---------------
Rt = TransposeVectorize(R) ; 
% Product of Rt and Xe. We first convert Rt in a block diagonal matrix 
Rt = ConvertBlockDiag(Rt) ; 
% 
Se = zeros(ndimB*nelemB,nnodeEb) ;
for inode=1:nnodeEb
    xeR = Rt*(Xe(:,inode)-Xe(:,1)) ;
    for idim = 1:ndimB
    Se(idim:ndimB:ndimB*nelemB,inode) = xeR(idim:ndim:ndim*nelemB) ;
    end
end