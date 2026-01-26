function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARCnew] = MAKE1ZERO_1STEP(xNEW,wNEW,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO)
%--------------------------------------------------------------------------
% function [xNEW, wNEW, CONVERGENCE, DATALOC, POLYINFO, VARCnew] = ...
%           MAKE1ZERO_1STEP(xNEW, wNEW, b, DATALOC, VAR_SMOOTH_FE, VARC, POLYINFO)
%
% PURPOSE:
%   Performs one Newton-based elimination step in the **first stage** of the 
%   Continuous Empirical Cubature Method (CECM).
%
%   It attempts to remove a single integration point by zeroing its weight,
%   and simultaneously updating the positions and weights of the remaining points 
%   to preserve the moment-matching condition:  Φᵗ·w ≈ b.
%
%   This step is repeated iteratively within `SPARSIF_1stStage` until the desired
%   number of cubature points is reached or until no further reduction is possible.
%
% INPUT:
%   - xNEW          : [n x ndim] Current coordinates of cubature points.
%   - wNEW          : [n x 1] Current weights.
%   - b             : [r x 1] Exact integrals of the basis functions: b = PHI' * W.
%   - DATALOC       : Struct with user-defined tolerances and control parameters:
%                     .TOL_NewtonRaphson_EliminationPoints
%                     .MaxIterationsNR_ElimPoints
%                     .MSGPRINT
%   - VAR_SMOOTH_FE : Struct with mesh, basis, and element data.
%   - VARC          : Struct defining active/inactive degrees of freedom:
%                     .POINTSl (free point DoFs), .POINTSRp (semi-fixed), etc.
%   - POLYINFO      : Struct with local polynomial coefficients, scaling, and triangulation.
%
% OUTPUT:
%   - xNEW          : Updated coordinates (possibly with one point removed).
%   - wNEW          : Updated weights (one weight set to zero).
%   - CONVERGENCE   : Flag (1: success; 0: unable to remove point).
%   - DATALOC       : Updated local data (including REMOVED_INDEX).
%   - POLYINFO      : Updated polynomial evaluation info (basis recomputed).
%   - VARCnew       : Updated control variable structure (reclassification of DoFs).
%
% METHOD OVERVIEW:
%   1. Select the candidate point to eliminate using a **significance criterion**:
%      S(i) = w_i * ||φ_i||².
%   2. Attempt to eliminate the selected point by applying constrained Newton iterations
%      (via `UPDATEPOSW`) to recompute the reduced cubature rule without it.
%   3. If convergence is successful (residual below tolerance), accept elimination.
%      Otherwise, proceed to the next candidate point.
%
% NOTES:
%   - The algorithm ensures that only positive weights remain at convergence.
%   - Points in `POINTSRp` are eventually relaxed into `POINTSl` after elimination.
%   - If no successful elimination is achieved, `CONVERGENCE = 0`.
%
% SEE ALSO:
%   - SPARSIF_1stStage
%   - UPDATEPOSW
%   - EVALBASIS
%   - EliminationStrategy
%
% REFERENCE:
%   J.A. Hernández et al., "Continuous Empirical Cubature Method for Data Compression
%   in Computational Mechanics", *Computer Methods in Applied Mechanics and Engineering*, 2024.
%
% AUTHOR:
%   Joaquín A. Hernández (UPC-CIMNE), April 2023.
%--------------------------------------------------------------------------

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp_48.mat')
end

DATALOC = DefaultField(DATALOC,'MSGPRINT',{}) ; % This is just for displaying information

[PHIk_y,~,POLYINFO]=     EVALBASIS(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO)  ;
errCOMP = norm(PHIk_y'*wNEW-b)./norm(b);      % MEASURE OF  fitting ERROR (when evaluation via fitting ) !!!
if DATALOC.iter == 1
    DATALOC.errorFITTING = errCOMP ;
    disp('--------------------------------------------') ;
    disp(['Integration appr. error at k= 0 -->',num2str(errCOMP*100),' %', '(  it measures the quality in evaluating the function via   fitting, if any)']);
end
% ------------------T
% Criterion for selecting the point to be removed: significance index
POINTS_F = [VARC.POINTSl(:); VARC.POINTSRp(:)] ; % All candidate points
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(POINTS_F)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL_NewtonRaphson_EliminationPoints ;
%-------------------------------
% Criterion for removing points
% ------------------------------
S_crit = wNEW(POINTS_F).*sum(PHIk_y(POINTS_F,:).*PHIk_y(POINTS_F,:),2) ;
[~, indREM ]= sort(S_crit) ;
xOLD = xNEW ;
wOLD = wNEW ;
VARCnew = VARC ;
% ----------------------------------------------------
iremove = 1 ;
SALIR = 0 ;
normB = norm(b);
%
while SALIR == 0 && iremove <=length(S_crit)
    % Elimination strategy
    % --------------------
    [VARCnew,xNEW,wNEW ]= EliminationStrategy(indREM,iremove,VARCnew,xNEW,wNEW) ;
    kMAX = DATALOC.MaxIterationsNR_ElimPoints  ;    k=1;    nF = 1e10 ;    CONVERGENCE = 1;
    nF_old = nF ;    nITER_increase = 0 ;    maxITER_increase = 4 ;
    ISOUT = 0 ;
    while  nF>=TOL && k<= kMAX && nITER_increase<=maxITER_increase
        [xkp1, wkp1, nF,POLYINFO,ISOUT,VARCnew  ] = ...
            UPDATEPOSW(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;
        nF = nF/normB ;    nd = (xkp1-xNEW).^2 ;    nd = max(sqrt(sum(nd,2))) ;
        disp(['Iteration k=',num2str(k),',  error residual =',num2str(nF),' MAX NORM incre DISPL =',num2str(nd),POLYINFO.MESSAGE_RANK])
        if ISOUT == 1
            break
        end
        DecreaseResidual = (nF_old-nF)/nF_old ;
        % Tolerance to consider that the residual is actually decreasing
        TOL_decrease_tolerance = 1e-6; % Set it to zero
        if  DecreaseResidual < TOL_decrease_tolerance
            nITER_increase = nITER_increase + 1 ;
        end
        xNEW = xkp1;
        wNEW = wkp1 ;
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
    if ( ISOUT ==1  ) || (k>kMAX  &&  nF>DATALOC.TOL_NewtonRaphson_EliminationPoints ) || nITER_increase >maxITER_increase
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;
        VARCnew = VARC ;
    else
        SALIR = 1 ;
        %   if  isempty(VARCnew.POINTSRpFIXED)
        % Points belonging to POINTSRp are moved to POINTSl
        VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)] ;
        VARCnew.POINTSRp = [] ;
        %         else
        %             % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables3D.mlx
        %             % Intersection between VARCnew.POINTSRp and VARCnew.POINTSRpFIXED
        %             VARCnew_POINTSRp_old = VARCnew.POINTSRp ;
        %             [VARCnew.POINTSRp,~,IndexesP] =   intersect(VARCnew.POINTSRpFIXED,VARCnew_POINTSRp_old) ;
        %             % Points which are not in the intersection are moved to VARCnew.POINTSl(:)
        %             IndexMoveToL = setdiff(1:length(VARCnew_POINTSRp_old),IndexesP) ;
        %             VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew_POINTSRp_old(IndexMoveToL)] ;
        %         end
    end
end
if  iremove >length(S_crit)
    CONVERGENCE = 0 ;
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ;
end
