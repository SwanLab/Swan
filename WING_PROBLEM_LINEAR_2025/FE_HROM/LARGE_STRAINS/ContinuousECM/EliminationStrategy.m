function [VARCnew,xNEW,wNEW ]= EliminationStrategy(indREM,iremove,VARCnew,xNEW,wNEW)
%--------------------------------------------------------------------------
% function [VARCnew, xNEW, wNEW] = EliminationStrategy(indREM, iremove, VARCnew, xNEW, wNEW)
%
% PURPOSE:
%   Implements the direct elimination strategy of a single integration point 
%   in the first-stage optimization of the Continuous Empirical Cubature Method (CECM).
%   The selected point has its weight set to zero, while its position remains fixed.
%
%   This step is part of a greedy iterative procedure for sparsifying 
%   the reduced integration rule while preserving the accuracy of the
%   integral approximation.
%
% INPUT:
%   - indREM    : [n x 1] Vector of candidate point indices (sorted by significance criterion).
%   - iremove   : Integer index into `indREM`, indicating which point to remove in this step.
%   - VARCnew   : Struct with the current classification of points:
%                .POINTSl   : Free points (weights & positions can change).
%                .POINTSRp  : Points with fixed positions, free weights.
%                .POINTSRpw : Points with fixed weights (=0) and positions (removed).
%   - xNEW      : [n x ndim] Coordinates of integration points.
%   - wNEW      : [n x 1] Weights associated to xNEW.
%
% OUTPUT:
%   - VARCnew   : Updated VARCnew structure reflecting the elimination.
%   - xNEW      : Unchanged coordinates of integration points.
%   - wNEW      : Updated weights with selected weight set to zero.
%
% FUNCTIONALITY:
%   - Moves the selected point from the set of free points (POINTSl)
%     to the set of fixed points with zero weight (POINTSRpw).
%   - Does **not** modify the point’s position (xNEW remains unchanged).
%
% REMARKS:
%   - This function assumes that `VARCnew.POINTSRp` is empty, as per the default
%     settings in the first-stage greedy elimination procedure.
%
% CONTEXT:
%   - Called from `MAKE1ZERO_1STEP` within the iterative loop of `SPARSIF_1stStage`.
%
% AUTHOR:
%   Joaquín A. Hernández (UPC-CIMNE), 2023.
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

%----------------------------------------------------

iremovLOC = indREM(iremove) ; % Indexes of point belonging to  POINTS_F which is constrained now (because we are going to set its weight to zero,
% and the position will remain unchanged)
POINTS_F = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)];
ipoint_control = POINTS_F(iremovLOC) ;  % Index of the constrained point (global)
% --------------------------------------------------------------------------------------
% We would have to guess now whether this point belongs to VARCnew.POINTSl or
% VARCnew.POINTSRp. However, in this strategy, VARCnew.POINTSRp is always
% empty. Therefore, we simply make 
VARCnew.POINTSl(iremovLOC)=[]  ;  % Remove the point from the unconstrained set
%POINTS_F=[VARCnew.POINTSl(:)'; VARCnew.POINTSRp(:)' ];  ;  % Remove the point from the unconstrained set

VARCnew.POINTSRpw = [VARCnew.POINTSRpw; ipoint_control] ;  % Add it to the weight/position contrained set 
% Now we must assign  to xNEW and wNEW the constrained values 
% In the case of xNEW, we do not need do to anything (it remains in the
% same position). 
% As for the weights, we simply make
wNEW(ipoint_control) = 0 ;  

disp('----------------------------------------')
disp(['Removing point =',num2str(ipoint_control),' (irem =',num2str(iremove),')'])
disp('---------------------------------------')