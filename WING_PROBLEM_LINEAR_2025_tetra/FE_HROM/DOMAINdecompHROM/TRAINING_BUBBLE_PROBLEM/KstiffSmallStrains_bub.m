function K = KstiffSmallStrains_bub(OPERFE,FgradST,ndim,celastST) 
if nargin == 0
    load('tmp.mat')
end
% Assembly Stiffness Matrix, large strains 
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 18
% 
% Compute celasLARGEgeo  
%celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ; 
% Compute celasLARGEmat 


if ~isempty(FgradST)
    celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;
    nF = size(celasLARGE,2) ;
    for icomp = 1:nF
        icol = icomp:nF:length(FgradST) ;
        celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ;
    end
    celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
    K = OPERFE.BstA'*(celasLARGE*OPERFE.BstA);
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
 %   celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;
 %
    nF = size(celastST,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
    end
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    K = OPERFE.BstA'*(celastST*OPERFE.BstA);
    
end
