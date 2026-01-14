function Bst_r = DetermineModeHomogViaREACT(DATAreactionsPLOT,MESH,Bst,W,OTHER_output)
% =========================================================================
% DETERMINEMODEHOMOGVIAREACT — Precompose reactions → homogenized PK1 map
% =========================================================================
% PURPOSE
%   Build a *precomposed* operator that lets the ROM go **directly from**
%   Gauss-point PK1 stresses to **homogenized PK1 components** (2D case).
%
%   Physically, the linear operator that turns boundary **reactions** into
%   homogenized stresses is **X_homog_stress**:
%
%       P_homog = X_homog_stress * R_bnd,
%
%   where R_bnd stacks the boundary nodal reactions (one DOF component per
%   boundary node). This routine forms the matrix product
%
%       Bst_r = Bst(:, DOFS_selected) * X_homog_stressᵀ,
%
%   which *precomposes* the extraction of boundary reaction DOFs (through
%   Bst and the selection DOFs_selected) with the geometric averaging
%   X_homog_stress. In other words, **X_homog_stress is the averaging map**,
%   while **Bst_r** adapts it to the ROM’s stress→residual assembly layout.
%
% WHAT IT DOES (2D implementation)
%   1) Collect all boundary nodes and their coordinates Xbnd = (X1, X2).
%   2) Build X_homog_stress ∈ ℝ^{4 × (2·N_bnd)} so that, for reactions
%      R_bnd = [R1₁,R2₁, R1₂,R2₂, …]ᵀ,
%          ⟨P11⟩ = Σ_i X1_i·R1_i ,  ⟨P22⟩ = Σ_i X2_i·R2_i
%          ⟨P12⟩ = Σ_i X2_i·R1_i ,  ⟨P21⟩ = Σ_i X1_i·R2_i
%      (optionally normalized by ΣW for volume/area averaging).
%   3) Expand boundary node ids → DOF ids: DOFS_selected = small2large(...).
%   4) Precompose with the ROM assembly matrix Bst:
%         Bst_r = Bst(:, DOFS_selected) * X_homog_stressᵀ.
%
%   Downstream, this allows you to use **Bst_r** in the same “stress-to-
%   forces” pipeline as other ROM/ECM operators, while keeping in mind that
%   **X_homog_stress is the true homogenization functional** acting on the
%   reaction vector.
%
% INPUTS
%   DATAreactionsPLOT : kept for interface consistency (unused here).
%   MESH :
%       .COOR        — nodal coordinates (to define X_homog_stress).
%       .NODES_FACES — cell with boundary faces (to collect boundary nodes).
%   Bst   : ROM assembly matrix consistent with internal-force construction
%           (stress → residual structure used elsewhere in the code).
%   W     : quadrature weights (used here for normalization by ΣW).
%   OTHER_output : reserved (e.g., for future normal-based formulations).
%
% OUTPUT
%   Bst_r : Precomposed operator compatible with the ROM assembly path.
%
% NOTES / LIMITATIONS
%   • Implemented for 2D (nSTRESS = 4: P11, P22, P12, P21). 3D extension TBD.
%   • Normalization by sum(W) is applied to X_homog_stress for averaging.
%   • Sign conventions correspond to PK1. Adapt X_homog_stress if you need
%     Cauchy/PK2 or different stress measures.
%
% VERSION / AUTHORSHIP
%   • 19-AUG-2025 — First version (Cartagena).
%   • 07-NOV-2025 — Header corrected: **X_homog_stress** is the linear
%                    homogenization operator; Bst_r is its precomposition
%                    with the ROM’s residual assembly. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

if nargin == 0
    load('tmp2.mat')
end


AllNodesBoundary = unique(cell2mat(MESH.NODES_FACES(:))) ; %{FACE_TO_ANALIZE} ;
Xbnd = MESH.COOR(AllNodesBoundary,:) ;

% Now we define a Mxnstress (nstress = 4 for  2D) containing
% rows that depends on Xbnd
if size(MESH.COOR,2) == 2
    nSTRESS = 4 ;
    X_homog_stress  = zeros(nSTRESS,2*size(Xbnd,1)) ;
    % PK1 homog. stress will be calculated as the sum of  size(Xbdn,1) matrices of the form
    % P_node_i = [React_1*X_1
    % React_2*X_2
    % React_1*X_2
    % React_2*X_1 ]
    %  = [X1   0
    %     0    X_2
    %     X2    0
    %     0     X1]    *
    % [React_1
    %  React_2]
    X_homog_stress(1,1:2:end) =    Xbnd(:,1) ;
    X_homog_stress(2,2:2:end) =    Xbnd(:,2) ;
    X_homog_stress(3,1:2:end) =    Xbnd(:,2) ;
    X_homog_stress(4,2:2:end) =    Xbnd(:,1) ;
    
else
    error('Option not implemented yet')
end

X_homog_stress = X_homog_stress/sum(W) ; 
 
ndim = size(MESH.COOR,2) ;
DOFS_selected = small2large(AllNodesBoundary,ndim) ;
 
% Resultant
Bst_r  = Bst(:,DOFS_selected)*X_homog_stress' ;
