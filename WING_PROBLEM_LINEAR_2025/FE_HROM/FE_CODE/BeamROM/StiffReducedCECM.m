function Kstiff_hyper = StiffReducedCECM(CECMoutput, nstrain,CgloDOM,DATA,BdomRED,Wdom,DATA_AUX)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/READMEmulti.mlx
if nargin == 0
    load('tmp.mat')
end

 

DATA = DefaultField(DATA,'ISNONLINEAR',0) ;

DATA_AUX = DefaultField(DATA_AUX,'Cglo_is_multiplied_by_the_weights',1) ; % = 0 ;
CelasBdomALL = [] ;

ngaus = DATA_AUX.MESH.ngausE ; 


if DATA.ISNONLINEAR == 0
    % HyperReduced elasticity matrix  ---including full-order weights
  
    CelasBdomALL = CgloDOM*BdomRED ;  % This plays the role of stresses
    for istrain = 1:nstrain       
        if DATA_AUX.Cglo_is_multiplied_by_the_weights == 1
           CelasBdomALL(istrain:nstrain:end,:) =  bsxfun(@times, CelasBdomALL(istrain:nstrain:end,:),1./Wdom) ;
        end
    end
    
    Celas_Bdom_hyper  =  InterpolationGaussVariablesECM(CelasBdomALL,CECMoutput,ngaus,nstrain) ;  
    BdomRED_hyper =  InterpolationGaussVariablesECM(BdomRED,CECMoutput,ngaus,nstrain) ;  
    wBdomRED_hyper = zeros(size(BdomRED_hyper)) ;
     for istrain = 1:nstrain
        INDS = istrain:nstrain:size(BdomRED_hyper,1) ;
        wBdomRED_hyper(INDS,:) = bsxfun(@times,BdomRED_hyper(INDS,:),CECMoutput.wCECM) ;
     end
    
    
    % Reduced stiffness matrix computed with reduced points
    Kstiff_hyper = wBdomRED_hyper'*Celas_Bdom_hyper  ;
    
else
    error('Option not implemented')
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
