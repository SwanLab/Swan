function [xNEW, wNEW, nF,POLYINFO,ISOUT,VARCnew  ] = UPDATEPOSW(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew)
%--------------------------------------------------------------------------
% function [xNEW, wNEW, nF, POLYINFO, ISOUT, VARCnew] = ...
%     UPDATEPOSW(wNEW, b, xNEW, DATALOC, VAR_SMOOTH_FE, POLYINFO, VARCnew)
%
% PURPOSE:
%   Performs one Newton-Raphson step to simultaneously update integration 
%   point positions (xNEW) and weights (wNEW) so as to reduce the error in 
%   approximating the integral of a reduced basis.
%
%   Specifically, solves the linearized system:
%       J * Δq ≈ b - PHI(xNEW)' * wNEW
%   where Δq stacks perturbations in positions and weights.
%
% INPUTS:
%   - wNEW           : [n x 1] Current integration weights.
%   - b              : [m x 1] Exact integral values (target vector).
%   - xNEW           : [n x ndim] Current coordinates of integration points.
%   - DATALOC        : Struct containing control parameters (tolerances, options).
%   - VAR_SMOOTH_FE  : Struct with mesh data, polynomial orders, etc.
%   - POLYINFO       : Struct containing cached basis function evaluations and metadata.
%   - VARCnew        : Struct with active DOF indices for weights and positions.
%
% OUTPUTS:
%   - xNEW           : Updated coordinates of integration points.
%   - wNEW           : Updated integration weights.
%   - nF             : Norm of the residual after update (measure of fitting error).
%   - POLYINFO       : Updated structure with basis info, rank messages, etc.
%   - ISOUT          : Flag = 1 if any point moved outside the domain mesh.
%   - VARCnew        : Updated DOF structure (e.g. DOFl).
%
% PROCEDURE:
%   1. Evaluates basis functions PHI and their gradients dPHI at xNEW.
%   2. Constructs the Jacobian J of the residual w.r.t. positions and weights.
%   3. Solves J * Δq = residual using rank-revealing SVD.
%   4. Applies update Δq and checks for domain exit.
%
% REMARKS:
%   - Only free DOFs (VARCnew.POINTSl) are included in the Jacobian.
%   - The solution Δq is dense unless J is full rank.
%   - The function assumes positions can be freely updated unless a point 
%     is part of the fixed DOF set (POINTSRp).
%
% AUTHOR:
%   Joaquín A. Hernández (UPC-CIMNE), 2023.
%--------------------------------------------------------------------------

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp.mat')
end
% Evaluation of the function and the gradients at all points%%
% --------------------------------------------------------
%setElementsBefore = POLYINFO.setElements; 
[PHIk_y,dPHIk_y,POLYINFO]=     EVALBASIS(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO);
% ------------------------------
bNEW = PHIk_y'*wNEW ;
ndim = size(xNEW,2) ;
% Residual
r = b-bNEW ;
nF = norm(r)  ;
% Computation of the Jacobian matrix (derivative with respect to X)
J = zeros(length(r),length(wNEW)*ndim);
for idim = 1:length(dPHIk_y)
    J(:,idim:ndim:end) =  bsxfun(@times,dPHIk_y{idim},wNEW)' ;
end
% We are only interested in the components associated to POINTSl
ndim = size(xNEW,2) ;
DOFl = small2large(VARCnew.POINTSl,ndim) ;
VARCnew.DOFl  = DOFl ;
J = J(:,DOFl) ;
% The other two blocks are formed by
POINTS_F = [VARCnew.POINTSl(:);VARCnew.POINTSRp(:) ] ;
J = [J,  PHIk_y(POINTS_F,:)'] ;
DATALOCSVD.RELATIVE_SVD = 1;
TOL = 1e-10 ;
[UU,SS,VV] = SVDT(J,TOL,DATALOCSVD) ;
% Solution of the underdetermined system of equations (more unknowns than equation)
if length(SS) == size(J,1)
    % No need to correct rank
    % -------------------------------------------------------
    delta_q = J\r ;    % Matlab produces a sparse solution (minimize  the l1 norm of delta_q)
    POLYINFO.MESSAGE_RANK = '' ;
else
   % disp(['Incomplete rank: number of equations = ',num2str(size(J,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']);
    UU = bsxfun(@times,UU',1./SS)'  ;
    delta_q = VV'\(UU'*r) ;
    POLYINFO.MESSAGE_RANK = ['(Rank ',num2str(length(SS)),' of ',num2str(size(J,1)),')'] ; 
end

[xNEW,wNEW,ISOUT,VARCnew,POLYINFO] =  UpdateAndCheckOutside2023(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;

