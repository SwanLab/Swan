function [M_hyper ]= MassReduced(setPoints, ngaus,densGLO,DATAIN,NstRED,WdomRED,Wdom,nstrain) ;

if nargin == 0
    load('tmp.mat')
    setPoints  =1:320*8 ; 
    WdomRED = Wdom ; 
end

% Computation  of the reduced stiffness using this set of points
% ---------------------------------------------------------------
%  Set of entries associated with these set of points
setDOFS = small2large(setPoints,nstrain) ;
% HyperReduced B matrix
NstRED_hyper = NstRED(setDOFS,:) ;
WdomRED = CompWeightDiag(WdomRED,nstrain) ;  ;

DATAIN = DefaultField(DATAIN,'ISNONLINEAR',0) ; 
if DATAIN.ISNONLINEAR == 0
    % HyperReduced elasticity matrix  ---including full-order weights
    densHYPER = densGLO(setDOFS,setDOFS) ;
    densHYPER = WdomRED*densHYPER ; 
    % Hyperreduced matrix multiplied by weights
  %  wNstRED_hyper = zeros(size(NstRED_hyper)) ;
    % PRoduct C*B
    M_hyper = (NstRED_hyper'*densHYPER)*NstRED_hyper ; 
    
%     Celas_Bdom = wCelas_hyper*NstRED_hyper ;
%     CelasBdom = zeros(size(wCelas_Bdom)) ;
% %     for istrain = 1:nstrain
% %         INDS = istrain:nstrain:length(setGauss) ;
% %         wNstRED_hyper(INDS,:) = bsxfun(@times,NstRED_hyper(INDS,:),WdomRED) ;
% %         Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
% %     end
%     
%     % Reduced stiffness matrix computed with reduced points
%     Kstiff_hyper = wNstRED_hyper'*Celas_Bdom  ;
    
else
    error('Option not implemented')
%     % NONLINEAR-PROBLEMS --- (J2-Plasticity )
%      % HyperReduced elasticity matrix  
%     Celas_hyper = densGLO(setGauss,setGauss) ;
%     % Hyperreduced matrix multiplied by weights
%     wNstRED_hyper = zeros(size(NstRED_hyper)) ;
%     % PRoduct C*B
%     Celas_Bdom = Celas_hyper*NstRED_hyper ;
%     
%     for istrain = 1:nstrain
%         INDS = istrain:nstrain:length(setGauss) ;
%         wNstRED_hyper(INDS,:) = bsxfun(@times,NstRED_hyper(INDS,:),WdomRED) ;
%       %  Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
%     end
%     
%     % Reduced stiffness matrix computed with reduced points
%     Kstiff_hyper = wNstRED_hyper'*Celas_Bdom  ;
    
    
end
