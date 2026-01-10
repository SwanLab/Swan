function Kc = KstiffLargeStrains_COROT_LRss(OPERFE,StwoST,FgradST,ndim,celastST,DATA,LboolCallQ,D_QrotALL,KcGEOunassemGLOloc)
%--------------------------------------------------------------------------
% function Kc = KstiffLargeStrains_COROT_LRss(OPERFE,StwoST,FgradST,ndim,...
%                                celastST,DATA,LboolCallQ,D_QrotALL,KcGEOunassemGLOloc)
%
% PURPOSE:
%   Computes the coarse-scale tangent stiffness matrix `Kc` for large strain
%   problems under a co-rotational framework. This implementation is adapted
%   to handle small strains combined with large rotations or large strains
%   with small rotations, within the framework of reduced-order models or
%   hyperreduction schemes (e.g., EIFEM).
%
%   The stiffness matrix is assembled as the sum of material and geometric
%   contributions, both evaluated in the current (rotated) configuration.
%
% INPUTS:
%   - OPERFE                : Structure containing precomputed reduced operators,
%                             weights, and matrices needed for assembling the stiffness.
%   - StwoST                : Second Piola–Kirchhoff stress tensor snapshots at CECM points.
%   - FgradST               : Deformation gradient snapshots at CECM points.
%   - ndim                  : Spatial dimension (2 or 3).
%   - celastST              : Material tangent moduli at the CECM points.
%   - DATA                  : Structure with global parameters and flags.
%   - LboolCallQ            : Boolean projection matrix (coarse-to-global dofs).
%   - D_QrotALL             : Rotation matrices applied at the coarse level.
%   - KcGEOunassemGLOloc    : Local unassembled geometric stiffness matrix.
%
% OUTPUT:
%   - Kc                    : Global coarse-scale tangent stiffness matrix,
%                             including material and geometric contributions.
%
% IMPLEMENTATION NOTES:
%   - The material component includes both the tangent stiffness due to material
%     nonlinearity and a geometric part derived from the second Piola–Kirchhoff
%     stress tensor.
%   - The geometric component (currently always included) represents internal
%     force effects induced by large displacements and rotations.
%   - The function handles both standard projection-based assembly and hyperreduced
%     forms, using local Gauss-point data and pre-assembled diagonal structures.
%
% REFERENCES:
%   - /home/joaquin/.../05_COROT_SSLR_LSSR.mlx
%   - EIFEM_largeROTfinal.pdf
%--------------------------------------------------------------------------

% Adaptation of KstiffLargeStrains.m to the co-rotational method
% See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 14-feb-2025, FRIDAY, 7:47, Balmes 185,  Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
if nargin == 0
    load('tmp1.mat')
end

% ---------------------------------------------------
% ------- MATERIAL COMPONENT, COARSE-SCALE ----------
% ---------------------------------------------------
% Compute celasLARGEgeo, FINE-SCALE, AT THE CECM POINTS
celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim) ;
% Compute celasLARGEmat, FINE-SCALE, AT THE CECM POINTS
celasLARGE = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;

% TOTAL
celasLARGE = celasLARGE + celasLARGEgeo ;



if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    error('Option not implemented')
    %     nF = size(celasLARGE,2) ;
    %
    %     for icomp = 1:nF
    %         icol = icomp:nF:length(FgradST) ;
    %         celasLARGE(icol,:) = bsxfun(@times,celasLARGE(icol,:),OPERFE.wSTs) ;
    %     end
    %
    %     celasLARGE = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
    %
    %     K = OPERFE.Bst'*(celasLARGE*OPERFE.Bst);
    
else
    %   \KcMAT = \LboolCallQ^T  \DiagC{\KcLINloc} \LboolCallQ +
    %  +  {\BmatIqST}^T \Par{\DiagC{\Wecm} (\DiagC{\CtangFst}-\DiagC{\CtangFlinST})} \BmatIqST
    %
    % where
    %  \BmatIqST = \DiagC{\BmatIst} \LboolCallQ
    BmatIqST = OPERFE.D_BmatIst*LboolCallQ ;
    % Linear component:  \LboolCallQ^T  \DiagC{\KcLINloc} \LboolCallQ
    Kc = LboolCallQ'*(OPERFE.D_KcLINloc*LboolCallQ);
    % Nonlinear component:        {\BmatIqST}^T \Par{\DiagC{\Wecm} (\DiagC{\CtangFst}-\DiagC{\CtangFlinST})} \BmatIqST
    Kc=  Kc +    BmatIqST'*(OPERFE.D_Wecm*((ConvertBlockDiag(celasLARGE)-OPERFE.D_CtangFlinST)*BmatIqST)) ;
    % ----------------------------------------------------------------------------------------------------------
    
end


% ---------------------------------------------------
% ------- GEOMETRIC COMPONENT, COARSE-SCALE ----------
% ---------------------------------------------------
% INCLUDE_GEO_MATRIX =  0;
% 
%  if INCLUDE_GEO_MATRIX == 1

Kc = Kc + OPERFE.LboolCall'*(KcGEOunassemGLOloc*LboolCallQ) ; 

 %end

% INCLUDE_GEO_MATRIX =  1; 
% 
% if INCLUDE_GEO_MATRIX == 1
% FintCunassemb  =  D_QrotALL*FintCunassembLOC ;
% D_FintCunassembGLO = DiagC_FintCunassemb(FintCunassemb,OPERFE.INDEXsparseFINT)  ;
% if DATA.MESH.ndim ==2
%     %     \KcGEO = \LboolCall^T   \DiagC{\AspMAT}  \DiagC{\FintCunassemb}   \DiagC{\PdownsRBlROT}   \LboolCallQ
%     KcGEO =  OPERFE.D_AspMAT*(D_FintCunassembGLO*OPERFE.D_PdownsRBlROT)  ; 
%     KcGEO = OPERFE.LboolCall'*(KcGEO*LboolCallQ)  ;
%     Kc = Kc +  KcGEO ; 
% else
%     error('Option not implemented')
% end
% end