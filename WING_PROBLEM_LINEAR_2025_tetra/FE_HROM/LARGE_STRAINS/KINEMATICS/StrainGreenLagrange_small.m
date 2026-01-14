function [Ev,FgradST ]= StrainGreenLagrange_small(gradU,ndim)
% Green-Lagrange strain tensor. Case Small strains

if nargin == 0
    load('tmp2.mat')
end


if ndim == 2
    % E = 0.5*(F'*F - I )
    nF = 4;
    nelem_ngaus = size(gradU,1)/nF  ;
    gradUcols_1 = 1:4:size(gradU,1) ;
    gradUcols_2 = 2:4:size(gradU,1) ;
    gradUcols_3 = 3:4:size(gradU,1) ;
    gradUcols_4 = 4:4:size(gradU,1) ;
    %
    nstrain = 3;
    Ev = zeros(nstrain*nelem_ngaus,size(gradU,2))  ;
    ECOLS_1 = 1:nstrain:size(Ev,1) ;
    ECOLS_2 = 2:nstrain:size(Ev,1) ;
    ECOLS_3 = 3:nstrain:size(Ev,1) ;
    
    Ev(ECOLS_1,:) = gradU(gradUcols_1,:);
    
    Ev(ECOLS_2,:) =  gradU(gradUcols_2,:) ;
    
    Ev(ECOLS_3,:) = gradU(gradUcols_3,:) +    gradU(gradUcols_4,:)   ;
    
else
    
    nF = 9 ;
    nelem_ngaus = length(gradU)/nF  ;
    gradUcols = cell(1,nF) ;
    for icols  =1:nF
        gradUcols{icols} = icols:nF:size(gradU,1) ;
    end
    nstrain = 6 ;
    Ev = zeros(nstrain*nelem_ngaus,size(gradU,2))  ;
    ECOLS = cell(1,nstrain) ;
    for icols = 1:nstrain
        ECOLS{icols} = icols:nstrain:size(Ev,1) ;
    end
    
    
    Ev(ECOLS{1},:) =  gradU(gradUcols{1},:) ;
    
    Ev(ECOLS{2},:) =  gradU(gradUcols{2},:)  ;
    
    Ev(ECOLS{3},:) = gradU(gradUcols{3},:)  ;
    
    Ev(ECOLS{4},:) =   gradU(gradUcols{4},:) +  gradU(gradUcols{7},:) ;
    
    Ev(ECOLS{5},:) =  gradU(gradUcols{5},:) + gradU(gradUcols{8},:)  ;
    
    Ev(ECOLS{6},:) =   gradU(gradUcols{6},:) +  gradU(gradUcols{9},:)  ;
end


% Now we set Fgrad = Identity
FgradST = zeros(size(gradU)) ;
for idim = 1:ndim
    LOCROWS = idim:ndim^2:size(FgradST,1) ;
    FgradST(LOCROWS,:) = 1;
end