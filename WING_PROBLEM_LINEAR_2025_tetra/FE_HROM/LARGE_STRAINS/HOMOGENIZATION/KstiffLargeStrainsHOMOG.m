function K = KstiffLargeStrainsHOMOG(OPERFE,StwoST,FgradST,ndim,celastST) ; 
% Assembly Stiffness Matrix, large strains 
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 18
% 
% Compute celasLARGEgeo  
celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ; 
% Compute celasLARGEmat 
celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ; 

celasLARGE = celasLARGE + celasLARGEgeo ; 

nF = size(celasLARGE,2) ; 

for icomp = 1:nF
    icol = icomp:nF:length(FgradST) ;
    celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ; 
end
 
celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix 

K = OPERFE.BstA'*(celasLARGE*OPERFE.BstA); 
