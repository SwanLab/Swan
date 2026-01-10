function [BstRED_r,DATAoffline] = DetermineModeReactionsOUTPUT_glo(DATAoffline,MESH,Bst,W,OTHER_output)
% =========================================================================
% DETERMINEMODEREACTIONSOUTPUT_GLO — Build reaction-resultant operator(s)
% =========================================================================
% PURPOSE
%   Provide the (optional) reduced operator BstRED_r that maps nodal
%   residuals to desired reaction RESULTANTS (JAHO: blunty speaking, it is not exactly like this ) to be considered in ECM/HROM.
%   We assume each target reaction resultant is a LINEAR COMBINATION of the
%   nodal residuals. This routine enables two mutually exclusive modes:
%
%   (A) FACE-BASED REACTIONS (local resultants on a selected boundary face)
%       • DATAoffline.DATAreactionsPLOT.FACE_TO_ANALIZE = [elemFaceList ...]
%       • Returns an operator BstRED_r that aggregates nodal residuals on
%         the requested face(s) into reaction components (Fx,Fy,Fz, Mx, …).
%
%   (B) HOMOGENIZED STRESSES VIA REACTIONS (global resultants)
%       • DATAoffline.DATAreactionsPLOT.ComputeHomogeneizedStresses = 1
%       • Returns BstRED_r that forms global “homogenized” stress-like
%         resultants from nodal residuals using quadrature weights W and
% >         mesh/BC info (OTHER_output).
%
%   If neither mode is active, or settings are inconsistent, the operator is
%   disabled (BstRED_r = []) and ECMforReactionForces is set to false.
%
% WHY THIS MATTERS
%   Including reaction resultants in A_fint (along with internal-force
%   snapshots) lets the ECM “see” boundary reaction content, improving
%   accuracy of hyperreduced models in problems where reactions matter
%   (e.g., constrained tests, homogenization, or reporting boundary forces).
%
% INPUTS
%   DATAoffline : Struct with reaction/ECM switches (fields set via DefaultField)
%       .ECMforReactionForces         (logical) master switch
%       .DATAreactionsPLOT            (struct or [])
%           .FACE_TO_ANALIZE              indices of boundary faces (mode A)
%           .ComputeHomogeneizedStresses  {0,1} homogenized mode (mode B)
%   MESH       : Mesh struct (connectivity, face sets, dimensions, etc.)
%   Bst        : Strain–displacement/assembly operator evaluated at Gauss pts
%                or the matrix that maps nodal unknowns to residual entries
%                used by the ROM (consistent with internal forces assembly).
%   W          : Quadrature weights vector/matrix for Gauss points (used in B)
%   OTHER_output : Extra FE/BC metadata (e.g., mapping matrices, sets, etc.)
%
% OUTPUTS
%   BstRED_r   : (matrix or []) Reaction-resultant operator to append to
%                BstRED_l*tau′(q) when building A_fint snapshots. Size is
%                compatible with stacking columns into A_fint.
%   DATAoffline: Possibly updated to ensure consistency:
%                   • .ECMforReactionForces = false when operator is not built.
%
% BEHAVIOR / LOGIC
%   • If ECMforReactionForces && DATAreactionsPLOT provided:
%       - If FACE_TO_ANALIZE is non-empty AND ComputeHomogeneizedStresses==0:
%           BstRED_r = DetermineModeReactionsOUTPUT(...)
%       - Else if ComputeHomogeneizedStresses==1:
%           (FACE_TO_ANALIZE ignored with a warning)
%           BstRED_r = DetermineModeHomogViaREACT(...)
%       - Else:
%           BstRED_r = []; ECMforReactionForces = false
%   • Else:
%       BstRED_r = []; ECMforReactionForces = false
%
% ASSUMPTIONS / NOTES
%   • Desired reactions are linear functionals of nodal residuals; this file
%     only selects/assembles the linear map. Actual residual content is built
%     elsewhere (stress → PK1 → internal forces).
%   • FACE_TO_ANALIZE should reference valid boundary faces; otherwise the
%     low-level builder will error.
%   • When homogenized stresses are requested, face selection is ignored and
%     a message is displayed.
%
% DEPENDENCIES (called internally)
%   DefaultField
%   DetermineModeReactionsOUTPUT         % builds face-based reaction resultants
%   DetermineModeHomogViaREACT           % builds homogenized-stress resultants
%
% VERSION / AUTHORSHIP
%   • 07-NOV-2025 — Comments added/clarified (global vs face-based modes, linear
%                    combination assumption, ECM integration). Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

if nargin == 0
    load('tmp1.mat')
end

DATAoffline = DefaultField(DATAoffline,'DATAreactionsPLOT',[]) ;
DATAoffline = DefaultField(DATAoffline,'ECMforReactionForces',0);
if DATAoffline.ECMforReactionForces  && ~isempty(DATAoffline.DATAreactionsPLOT)
    
    DATAoffline.DATAreactionsPLOT = DefaultField(DATAoffline.DATAreactionsPLOT,'ComputeHomogeneizedStresses',0) ;
    DATAoffline.DATAreactionsPLOT = DefaultField(DATAoffline.DATAreactionsPLOT,'FACE_TO_ANALIZE',[]) ;
    if ~isempty(DATAoffline.DATAreactionsPLOT.FACE_TO_ANALIZE) && DATAoffline.DATAreactionsPLOT.ComputeHomogeneizedStresses == 1
        disp(['Both ComputeHomogeneizedStresses and FACE_TO_ANALIZE are active, yet we will only compute  HOMOGENEIZED STRESSES'])
        DATAoffline.DATAreactionsPLOT.FACE_TO_ANALIZE = [] ;
    end
    if ~isempty(DATAoffline.DATAreactionsPLOT.FACE_TO_ANALIZE )
        BstRED_r = DetermineModeReactionsOUTPUT(DATAoffline.DATAreactionsPLOT,MESH,Bst);
        
    else
        if DATAoffline.DATAreactionsPLOT.ComputeHomogeneizedStresses == 1
            BstRED_r = DetermineModeHomogViaREACT(DATAoffline.DATAreactionsPLOT,MESH,Bst,W,OTHER_output);
        else
            BstRED_r = [] ;
            DATAoffline.ECMforReactionForces = false;
        end
    end
    
else
    BstRED_r = [] ;
    DATAoffline.ECMforReactionForces  =false;
end