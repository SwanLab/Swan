function Se = ChangeCoordBndVect(Xe,ndimB)
% Given the   matrix (Xe) of global coordinates of the position of nodes of boundary element "e", this function
% expresses these positions in a coordinate system intrinsic to the
% element (Se).
% Vectorized format
% ------------------
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 27-Oct-2015
if nargin == 0
    load('tmp.mat')
end
nnodeEb = size(Xe,2) ;
ndim = ndimB+1 ;
nelemB = size(Xe,1)/ndim ;  


% Rotation matrix 
R = zeros(nelemB*ndim,ndim) ; 
% Vector t1
t1 = Xe(:,2)-Xe(:,1) ;   
t1 = NormalizeVector(t1,ndim) ; % Normalization
if ndimB==1
    ind2 = 2:ndim:ndim*nelemB ; 
    ind1 = 1:ndim:ndim*nelemB ; 
    R(:,1) = t1 ;   % ERROR DETECTED IN PREVIOUS VERSIONS OF THE CODE. tHIS NEW VERSION IS OK
    R(ind1,2)=[-t1(ind2) ];  
     R(ind2,2)=[ t1(ind1)];
elseif ndimB == 2
    R(:,1) = t1 ;
    t2 = Xe(:,3)-Xe(:,1) ;    t2 =NormalizeVector(t2,ndim) ;
    %% Cross-product 
     t1t2 = CrossProductVect(t1,t2,ndim) ;     
    R(:,3) =  NormalizeVector(t1t2,ndim)  ; 
    if any(isnan(R(:,3)))
        error('Some of the chosen points are aligned')
    end
    R(:,2) = CrossProductVect(R(:,3),R(:,1),ndim) ;    
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