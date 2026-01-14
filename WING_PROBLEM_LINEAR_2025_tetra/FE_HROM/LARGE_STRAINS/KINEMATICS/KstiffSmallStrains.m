function K = KstiffSmallStrains(OPERFE,FgradST,ndim,celastST,DATA)
%--------------------------------------------------------------------------
% K = KstiffSmallStrains(OPERFE, FgradST, ndim, celastST, DATA)
%
% PURPOSE:
%   Assemble the global stiffness matrix under small strain assumptions,
%   using either the consistent tangent modulus or the linearized version
%   adapted for Empirical Cubature or Empirical Interpolation contexts (ECM/CECM).
%
% DESCRIPTION:
%   Depending on the availability of the deformation gradient (FgradST), and
%   the value of DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES, this routine computes:
%     - The full consistent stiffness matrix using Bᵀ·C·B formulation, with C 
%       either as material or geometrically nonlinear tangent.
%     - An ECM-based modified stiffness matrix, which corrects the tangent 
%       operator using an offline-stored reference configuration.
%
%   The tangent operator C is built using `CelasLARGEmat_allgauss`, and optionally
%   postprocessed using tailored weights or ECM corrections.
%
% INPUT:
%   OPERFE     : Structure containing finite element operators:
%                - Bst  : strain-displacement matrix at Gauss points
%                - wSTs : Gauss weights
%                - celastST_ini (optional): reference tangent operator (CECM)
%                - KstiffLINEAR (optional): linear part for CECM correction
%
%   FgradST    : Deformation gradient at Gauss points (can be empty for linear regime)
%   ndim       : Spatial dimension (2 or 3)
%   celastST   : Elastic tangent operator at Gauss points (Voigt format)
%   DATA       : Structure containing solver settings, including:
%                - CECM_ONLY_FOR_NONLINEAR_STRESSES : flag for inelastic correction
%
% OUTPUT:
%   K          : Global stiffness matrix assembled with weighted tangent operator
%
% NOTES:
%   - When `FgradST` is empty, stiffness is built from `celastST` only.
%   - If `DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1`, then the tangent
%     is corrected using stored initial tangents (`OPERFE.celastST_ini`).
%   - Tailored weights via `OPERFE.BstW` are no longer used (commented out).
%
% RELATED FILES:
%   - CelasLARGEmat_allgauss.m : Material tangent operator assembly
%   - ConvertBlockDiag.m       : Converts blockwise Gauss-wise tangent to matrix
%   - README_RigidBodyMotions.pdf (Page 18)
%   - 19_ExactLinearStiff.mlx (for CECM correction)
%--------------------------------------------------------------------------

% Assembly Stiffness Matrix, large strains
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 18
%
% Compute celasLARGEgeo
%celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ;
% Compute celasLARGEmat

if nargin == 0
    load('tmp1.mat')
end


if ~isempty(FgradST)
    
    if ~isfield(OPERFE,'BstW')
        
        celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;
        nF = size(celasLARGE,2) ;
        for icomp = 1:nF
            icol = icomp:nF:length(FgradST) ;
            celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ;
        end
        celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
        K = OPERFE.Bst'*(celasLARGE*OPERFE.Bst);
        
    else
        
        % Tailored weights for the ECM, see
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
        % This option turned out to be unreliable ! 
        celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;
        
        celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
        K = OPERFE.BstW'*(celasLARGE*OPERFE.Bst);
        
    end
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
    %   celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;
    %
    
  %  if ~isfield(OPERFE,'BstW')
        
       
   if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0              
              
        
        nF = size(celastST,2) ;
        for icomp = 1:nF
            icol = icomp:nF:size(celastST,1) ;
            celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
        end
        celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
        K = OPERFE.Bst'*(celastST*OPERFE.Bst);
        
   else       
       %------------------------------------------------------------------------------------------------------------------------
       % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
       % EIFEM, treating separaterly inelastic stresses, see also 
       % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/InternalForces.m
       % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/CONSTITUTIVE_MODELS/PK2stress_Constitutive_ModelVAR.m
       %-------------------------------------------------------------------------------------------------------------------------
         nF = size(celastST,2) ;
        for icomp = 1:nF
            icol = icomp:nF:size(celastST,1) ;
            celastST(icol,:) = bsxfun(@times,(celastST(icol,:)-OPERFE.celastST_ini(icol,:)),OPERFE.wSTs) ;
        end
        celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
        K = OPERFE.KstiffLINEAR + OPERFE.Bst'*(celastST*OPERFE.Bst);
       
   end
        
        
   % else        
          % Tailored weights for the ECM, see
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
    % OPTION NO LONGER USED ! YOU MAY DELETE IT
    %    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
     %   K = OPERFE.BstW'*(celastST*OPERFE.Bst);
        
   % end
    
end
