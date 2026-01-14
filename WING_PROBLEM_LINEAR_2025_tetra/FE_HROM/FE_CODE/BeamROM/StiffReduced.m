function [Kstiff_hyper,Celas_Bdom,CelasBdomALL ]= StiffReduced(setPoints, nstrain,CgloDOM,DATA,BdomRED,WdomRED,Wdom,DATA_GENGAUSS)

if nargin == 0
    load('tmp1.mat')
end

% Computation  of the reduced stiffness using this set of points
% ---------------------------------------------------------------
%  Set of entries associated with these set of points
setGauss = small2large(setPoints,nstrain) ;
% HyperReduced B matrix
BdomRED_hyper = BdomRED(setGauss,:) ;

DATA = DefaultField(DATA,'ISNONLINEAR',0) ;

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'Cglo_is_multiplied_by_the_weights',1) ; % = 0 ;
CelasBdomALL = [] ;


if DATA.ISNONLINEAR == 0
    % HyperReduced elasticity matrix  ---including full-order weights
    wCelas_hyper = CgloDOM(setGauss,setGauss) ;
    % Hyperreduced matrix multiplied by weights
    wBdomRED_hyper = zeros(size(BdomRED_hyper)) ;
    % PRoduct C*B
    wCelas_Bdom = wCelas_hyper*BdomRED_hyper ;
    Celas_Bdom = wCelas_Bdom ;
    CelasBdomALL = CgloDOM*BdomRED ;  % For reconstruction purposes
    for istrain = 1:nstrain
        INDS = istrain:nstrain:length(setGauss) ;
        wBdomRED_hyper(INDS,:) = bsxfun(@times,BdomRED_hyper(INDS,:),WdomRED) ;
        if DATA_GENGAUSS.Cglo_is_multiplied_by_the_weights == 1
            Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
            CelasBdomALL(istrain:nstrain:end,:) =  bsxfun(@times, CelasBdomALL(istrain:nstrain:end,:),1./Wdom) ;
        end
    end
    
    % Reduced stiffness matrix computed with reduced points
    Kstiff_hyper = wBdomRED_hyper'*Celas_Bdom  ;
    
else
    
    % NONLINEAR-PROBLEMS --- (J2-Plasticity )
    % HyperReduced elasticity matrix
    Celas_hyper = CgloDOM(setGauss,setGauss) ;
    % Hyperreduced matrix multiplied by weights
    wBdomRED_hyper = zeros(size(BdomRED_hyper)) ;
    % PRoduct C*B
    Celas_Bdom = Celas_hyper*BdomRED_hyper ;
    
    for istrain = 1:nstrain
        INDS = istrain:nstrain:length(setGauss) ;
        wBdomRED_hyper(INDS,:) = bsxfun(@times,BdomRED_hyper(INDS,:),WdomRED) ;
        %  Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
    end
    
    % Reduced stiffness matrix computed with reduced points
    Kstiff_hyper = wBdomRED_hyper'*Celas_Bdom  ;
    
    
end
