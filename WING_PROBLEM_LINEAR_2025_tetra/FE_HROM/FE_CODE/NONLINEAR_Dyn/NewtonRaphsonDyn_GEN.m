function [GAUSSV_np1,NODESV_np1,iter,CONVERGENCE] = ...
    NewtonRaphsonDyn_GEN(DATA,NODESV_n,GAUSSV_n,BBfNW,...
    istep,F_f,OPERfe,B_s_g0,...
    PROPMAT,massMf,massMs,BND_Dyn,DynDATA,increT,BBf,wST,ASSEMBLY_INFO) ;

if nargin == 0
    load('tmp1.mat')
end
% GEneral function for Newton-Raphson
% ------------------------- JAHO, 24-Sept.2018

%%
DATA = DefaultField(DATA,'TypeImplementation','J2_plasticity_small_strains') ;
DATA  = DefaultField(DATA,'tolNEWTONRAPSHON',1e-6) ;
DATA = DefaultField(DATA,'max_iter',20) ;



switch  DATA.TypeImplementation
    case 'J2_plasticity_small_strains'
        % Only function implemented so far, for 2D (24-Sept-2018)
        % This function must be made more general in future implementations
        % ---------------------------
        if ~isempty(BND_Dyn.gBOUNDdd)
        gBOUNDdd = BND_Dyn.gBOUNDdd(:,istep) ; 
        else
            gBOUNDdd = 0 ; 
        end
        [GAUSSV_np1.stressST,NODESV_np1.U,...
            NODESV_np1.VEL(OPERfe.DOFf),NODESV_np1.ACEL(OPERfe.DOFf),...   %  NODESV_np1.VEL(OPERfe.DOFf),NODESV_np1.ACEL(OPERfe.DOFf),...
            GAUSSV_np1.EP,GAUSSV_np1.alpha,...
            GAUSSV_np1.sigmay,iter, CONVERGENCE,TIMEelap] =...
            NewtonRaphsonDyn_FE2D(DATA.tolNEWTONRAPSHON,...
            GAUSSV_n.EP,GAUSSV_n.sigmay,GAUSSV_n.alpha,...
            BBfNW,istep ,NODESV_n.U,DATA.max_iter,F_f,OPERfe.DOFf,B_s_g0,...
            PROPMAT,  ...
            massMf,massMs,gBOUNDdd,DynDATA,NODESV_n.VEL(OPERfe.DOFf),...
            NODESV_n.ACEL(OPERfe.DOFf),increT,BBf',wST,ASSEMBLY_INFO) ;
        
        % Compute Von Mises 
        % --------------------
        [ GAUSSV_np1.VonMises ] =  VonMises_Stress(reshape(GAUSSV_np1.stressST,OPERfe.nstrain,[] ))' ;
        
        
end

