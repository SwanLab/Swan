function [BasisF,SingVal_intf,nstrain,SNAPforceSnw,BasisS,V,DATAOUT] = ...
    BasisFfromStress(BasisS,SingVal_stress,BdomRED, Wdom,DATAIN)
% See Implementation.pdf, section MAtrix of internal forces
%dbstop('4')
if nargin ==0
    load('tmp.mat')
    % DATACUB = [] ;
end
DATAIN = DefaultField(DATAIN,'DATACUB',[]) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
%DATAIN =    DefaultField(DATAIN,'NOTMULTIPLIED_BY_WEIGHTS',0) ;
DATACUB = DATAIN.CUBATURE ;

DATAIN = DefaultField(DATAIN,'BasisS_is_directly_snapshotsINTERNALFORCES',0) ;

if DATAIN.BasisS_is_directly_snapshotsINTERNALFORCES == 0
    DATACUB= DefaultField(DATACUB,'IncludeSingularValuesStress',1) ; % Include Singular values stresses
    if DATACUB.IncludeSingularValuesStress == 1
        BasisS = bsxfun(@times,BasisS',SingVal_stress)' ;
    end
end
ngaus = length(Wdom);
nstrain = size(BdomRED,1)/length(Wdom);

DATAIN = DefaultField(DATAIN,'IMPOSE_VOLUME_CONSTRAINT',0) ;

% Product BasisRED^T*BasisS

if DATAIN.BasisS_is_directly_snapshotsINTERNALFORCES == 0
    % --------------------------------------------------------
    SNAPforceS = BasisRED_BasisStress(BdomRED,BasisS,nstrain,ngaus) ;
    %--------------------------------------------------------------
else
    SNAPforceS = BasisS ;
end
SNAPforceSnw = SNAPforceS ;

DATAOUT = [] ;

if ~exist('RSVDT')
    addpath('SVDlibrary')
end



if DATAIN.IMPOSE_VOLUME_CONSTRAINT ==1
    % Columns
    SNAPforceSnw = SNAPforceS ;
    intEXACT = SNAPforceS'*Wdom ;
    VOL = sum(Wdom) ;
    SNAPforceS = bsxfun(@minus,SNAPforceS',intEXACT/VOL)' ;
    SNAPforceS = bsxfun(@times,SNAPforceS,sqrt(Wdom)) ;
    
elseif DATAIN.IMPOSE_VOLUME_CONSTRAINT ==0  ||  DATAIN.IMPOSE_VOLUME_CONSTRAINT ==3
    SNAPforceS = bsxfun(@times,SNAPforceS,sqrt(Wdom)) ;
elseif  DATAIN.IMPOSE_VOLUME_CONSTRAINT ==2
    % We add another column
    FactorFirstColumnSNAPforce  =norm(SNAPforceS,'fro')/size(SNAPforceS,2) ;
    % including the square of the weights (multiplied by nnn), to make
    % the column "dominant" in an statistical sense
    % In doing so, we ensure tha the sum of the weights is equal to
    % the volume
    SNAPforceS = [FactorFirstColumnSNAPforce*ones(size(Wdom)),  SNAPforceS] ;
    SNAPforceS = bsxfun(@times,SNAPforceS,sqrt(Wdom)) ;
    DATAOUT.FactorFirstColumnSNAPforce = FactorFirstColumnSNAPforce ;
else
    error('Option not implemented')
end




hold on
nfigure = 100;
LEGENDG = 'Internal Virtual Work ' ;
COLOR = 'k-' ;
%dbstop('40')
DATAIN = DefaultField(DATAIN,'TOL_LOC_InternalForces',1e-6) ;
DATAIN = DefaultField(DATAIN,'TOL_LOC_InternalForces_IS_INTEGRATION',0) ;  % 9th-April-2020 (27th day quarantine COVID19)


if DATAIN.TOL_LOC_InternalForces_IS_INTEGRATION == 0
    % TRUNCATION CRITERION DICTATED BY FROBENIOUS NORM (STANDARD SVD)
    DATAIN.TOL_LOC = DATAIN.TOL_LOC_InternalForces  ;
    [U,S,V,h3,h4] = SVD_and_error(SNAPforceS,nfigure,LEGENDG,[],COLOR,DATAIN ) ;
    hold on
    legend([h3 ],{'Int. virt. work'})
    legend([h4 ],{'Int. virt. work'})
    BasisF = U  ; % = [2] ;
    SingVal_intf = S ;
    
    %     SNAPforceS_approx = bsxfun(@times,U',S)'*V';
    %     E= SNAPforceS-SNAPforceS_approx ;
    %     save('tmp.mat','E')
    
else
    
    
    
    % Truncation criterion based on integration error
    [U,S,V] = RSVDT(SNAPforceS) ;
    intEXACT = sqrt(Wdom)'*SNAPforceS ;
    [BasisF] = BasisF_truncationINTEGRATION(intEXACT,U,S,V,Wdom,DATAIN) ;
    SingVal_intf  =ones(size(BasisF,2),1) ;
    
end


if DATAIN.IMPOSE_VOLUME_CONSTRAINT ==1
    BasisF = [sqrt(Wdom),BasisF]  ;
    SingVal_intf = ones(size(BasisF,2),1) ;    % Not important, just for monitoring the error in the ECM
elseif DATAIN.IMPOSE_VOLUME_CONSTRAINT ==3
    % This option is appropriate when the integrals of the modes are all
    % zero
    BasisF = [sqrt(Wdom),BasisF]  ;
    
end





