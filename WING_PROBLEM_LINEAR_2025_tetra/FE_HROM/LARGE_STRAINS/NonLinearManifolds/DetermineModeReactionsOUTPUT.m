function Bst_r = DetermineModeReactionsOUTPUT(DATAreactionsPLOT,MESH,Bst)
% =========================================================================
% DETERMINEMODEREACTIONSOUTPUT — Face-based reaction resultant operator
% =========================================================================
% PURPOSE
%   Build a LINEAR FUNCTIONAL that maps nodal residuals on a selected mesh
%   FACE to a scalar REACTION RESULTANT. We assume the target reaction is a
%   linear combination of the nodal residual entries (e.g., sum of forces
%   along a chosen DOF on that face).
%
% HOW IT WORKS
%   1) Pick the boundary face: DATAreactionsPLOT.FACE_TO_ANALIZE.
%   2) Gather nodes on that face and expand to DOFs via small2large.
%   3) Select a single component per node using DATAreactionsPLOT.DOF_TO_ANALIZE
%      (1=x, 2=y, 3=z for 3D).
%   4) Form the linear combiner by summing those residual DOF rows:
%        Bst_r = Bst(:, DOFS_selected) * ones(#selected,1)
%      so, when applied to nodal unknowns/residuals (ACTUALLY STRESSES!), it produces the desired
%      scalar reaction resultant.
%
% INPUTS
%   DATAreactionsPLOT :
%       .FACE_TO_ANALIZE   — index of the boundary face to analyze.
%       .DOF_TO_ANALIZE    — component to aggregate per node (1..ndim).
%   MESH :
%       .NODES_FACES{f}    — node list for face f.
%       .COOR              — coordinates (used to infer ndim).
%   Bst  : Assembly operator consistent with the nodal residual vector used
%          in the ROM (same ordering as the global DOF vector).
%
% OUTPUTS
%   Bst_r : [nres × 1] operator that linearly aggregates nodal residual DOFs
%           on the selected face into a scalar reaction resultant.
%
% ASSUMPTIONS / NOTES
%   • Resultant is a plain sum. If you need face-measure weighting, signs,
%     or traction orientation, incorporate them in Bst or replace the final
%     ones(…) with appropriate weights.
%   • For vector reactions (e.g., Fx and Fy), call this routine per DOF
%     component and horizontally concatenate the resulting columns.
%   • FACE_TO_ANALIZE must be valid and reference a boundary face whose node
%     list is in MESH.NODES_FACES.
%
% DEPENDENCIES
%   small2large  — expands node indices to DOF indices (per spatial dim).
%
% VERSION / AUTHORSHIP
%   • 07-NOV-2025 — Initial documented header for face-based scalar reaction
%                    operator. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================
 

FACE_TO_ANALIZE = DATAreactionsPLOT.FACE_TO_ANALIZE ;
NODES_FACES= MESH.NODES_FACES{FACE_TO_ANALIZE} ;
ndim = size(MESH.COOR,2) ;
DOFs_face = small2large(NODES_FACES,ndim) ;
DOFS_selected = DOFs_face(DATAreactionsPLOT.DOF_TO_ANALIZE:ndim:end) ;

% Resultant 
Bst_r  = Bst(:,DOFS_selected)*ones(length(DOFS_selected),1) ; 
 