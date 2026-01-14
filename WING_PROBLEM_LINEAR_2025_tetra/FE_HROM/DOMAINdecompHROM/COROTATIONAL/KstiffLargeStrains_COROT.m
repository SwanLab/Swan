function K = KstiffLargeStrains_COROT(OPERFE,StwoST,FgradST,ndim,celastST,DATA,KcLINq,BstQ)
% Adaptation of KstiffLargeStrains.m to the co-rotational method
% JAHO, 29-Oct-2024, Tuesday, Honest Greens, Pedralbes Center, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
if nargin == 0
    load('tmp1.mat')
end
 
%
% Compute celasLARGEgeo
celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ;
% Compute celasLARGEmat
celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;

celasLARGE = celasLARGE + celasLARGEgeo ;



if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    error('Option not implemented')
    nF = size(celasLARGE,2) ;
    
    for icomp = 1:nF
        icol = icomp:nF:length(FgradST) ;
        celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ;
    end
    
    celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
    
    K = OPERFE.Bst'*(celasLARGE*OPERFE.Bst);
    
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
%     nF = size(celasLARGE,2) ;
%     
%     for icomp = 1:nF
%         icol = icomp:nF:length(FgradST) ;
%         celasLARGE(icol,:) = bsxfun(@times,(celasLARGE(icol,:)-OPERFE.celastST_ini(icol,:)),OPERFE.wSTs) ;
%     end
%     
%     celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
%     
%     K =OPERFE.KstiffLINEAR +  OPERFE.Bst'*(celasLARGE*OPERFE.Bst);

% NEW APPRROACH, EIFEM
%   \Kc =   \KcLINq + \BstQ^T \DiagC{\Wecm} ( \DiagC{\CtangFst} - \DiagC{\CtangFlinST})   \BstQ 
% ----------------------------------------------------------------------------------------------
K=   BstQ'*(OPERFE.D_Wecm*((ConvertBlockDiag(celasLARGE)-OPERFE.D_CtangFlinST)*BstQ)) ; 
K = K + KcLINq ; 
    
end

