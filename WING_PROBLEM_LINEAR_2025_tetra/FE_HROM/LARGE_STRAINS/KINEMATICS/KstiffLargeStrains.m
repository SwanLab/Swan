function K = KstiffLargeStrains(OPERFE,StwoST,FgradST,ndim,celastST,DATA)
if nargin == 0
    load('tmp.mat')
end
% ---------------------------------------------------------------------------------------------------
% FUNCTION: KstiffLargeStrains
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Assembles the **global tangent stiffness matrix K** for a hyperelastic problem under **large strains**
%   using a total Lagrangian formulation. It combines the **material tangent moduli** and the **geometric stiffness**
%   (stress-dependent) contributions across all Gauss points.
%
%   The final matrix takes the form:
%       K = K_mat + K_geo
%   where:
%       - K_mat accounts for nonlinear material behavior (via `celastST`)
%       - K_geo accounts for stress-induced stiffness (via `StwoST`)
%
% USAGE:
%   K = KstiffLargeStrains(OPERFE, StwoST, FgradST, ndim, celastST, DATA)
%
% INPUTS:
%   - OPERFE     : Structure with finite element operators, including:
%                   * Bst: strain-displacement matrices
%                   * wSTs: Gauss weights
%                   * celastST_ini: initial tangent tensor (optional)
%                   * KstiffLINEAR: reference stiffness (if linear part precomputed)
%   - StwoST     : Second Piola–Kirchhoff stress tensor at Gauss points (vectorized)
%   - FgradST    : Deformation gradient tensor at Gauss points (vectorized)
%   - ndim       : Problem dimension (2D or 3D)
%   - celastST   : Material tangent stiffness tensor (vectorized)
%   - DATA       : Structure with solver settings (e.g., `CECM_ONLY_FOR_NONLINEAR_STRESSES` flag)
%
% OUTPUT:
%   - K : Global stiffness matrix assembled from Gauss-point contributions
%
% FUNCTIONALITY:
%   - Computes geometric stiffness matrix using `CelasLARGEgeo_allgauss`.
%   - Computes material tangent stiffness (consistent tangent operator) using `CelasLARGEmat_allgauss`.
%   - If `DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0`:
%       * Fully assembles total stiffness from scratch.
%   - If `DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1`:
%       * Uses precomputed linear stiffness and adds only nonlinear part via incremental `celasLARGE - celastST_ini`.
%   - Applies Gauss quadrature weights to each block row of the tangent matrix.
%   - Uses `ConvertBlockDiag` to assemble block-diagonal format from Gauss-wise tangent blocks.
%
% REFERENCES:
%   - Derivation and formulation:
%     /LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/README_RigidBodyMotions.pdf, page 18
%   - Use case in:
%     /TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE  
%   Date: 6-Jan-2024 (adapted for CECM nonlinearity separation)
%   Comments by ChatGPT4, 12-May-2025
% DEPENDENCIES:
%   - CelasLARGEmat_allgauss
%   - CelasLARGEgeo_allgauss
%   - ConvertBlockDiag
%   - bsxfun (for applying Gauss weights efficiently)
%
% ---------------------------------------------------------------------------------------------------

% Assembly Stiffness Matrix, large strains
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 18
%
% Compute celasLARGEgeo
celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ;
% Compute celasLARGEmat
celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;

celasLARGE = celasLARGE + celasLARGEgeo ;



if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    
    nF = size(celasLARGE,2) ;
    
    for icomp = 1:nF
        icol = icomp:nF:length(FgradST) ;
        celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ;
    end
    
    celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
    
    K = OPERFE.Bst'*(celasLARGE*OPERFE.Bst);
    
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
    nF = size(celasLARGE,2) ;
    
    for icomp = 1:nF
        icol = icomp:nF:length(FgradST) ;
        celasLARGE(icol,:) = bsxfun(@times,(celasLARGE(icol,:)-OPERFE.celastST_ini(icol,:)),OPERFE.wSTs) ;
    end
    
    celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
    
    K =OPERFE.KstiffLINEAR +  OPERFE.Bst'*(celasLARGE*OPERFE.Bst);
    
end

