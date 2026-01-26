function  [xOLD,wOLD,DATALOC,POLYINFO,VARC,HISTORY] =...
    SPARSIF(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,POLYINFO)
%--------------------------------------------------------------------------
% function [xOLD, wOLD, DATALOC, POLYINFO, VARC, HISTORY] = ...
%     SPARSIF(xOLD, wOLD, b, DATALOC, VAR_SMOOTH_FE, VARC, HISTORY, POLYINFO)
%
% PURPOSE:
%   Performs the second stage of the Continuous Empirical Cubature Method (CECM),
%   applying a progressive, iterative elimination of integration points.
%   In each iteration, one point is removed and the weights/positions are 
%   updated via a constrained optimization step to preserve the accuracy 
%   of the quadrature rule.
%
% DESCRIPTION:
%   This function progressively sparsifies the integration rule computed 
%   during the first stage (direct elimination), further reducing the 
%   number of integration points while controlling the quadrature error.
%
% INPUTS:
%   - xOLD           : [n x ndim] Initial coordinates of integration points.
%   - wOLD           : [n x 1] Initial integration weights.
%   - b              : [m x 1] Exact integral of the basis functions.
%   - DATALOC        : Struct with control parameters (tolerances, limits, etc.).
%   - VAR_SMOOTH_FE  : Struct containing mesh and polynomial degree info.
%   - VARC           : Struct with DOF classification for weights and positions.
%   - HISTORY        : Struct storing the evolution of points/weights over iterations.
%   - POLYINFO       : Struct with metadata for basis function evaluation.
%
% OUTPUTS:
%   - xOLD           : Final coordinates of the selected ECM points.
%   - wOLD           : Final optimized weights.
%   - DATALOC        : Updated control parameters (with iteration data).
%   - POLYINFO       : Updated basis function metadata.
%   - VARC           : Updated DOF structure.
%   - HISTORY        : Updated with second-stage iteration history.
%
% PROCEDURE:
%   1. Initializes the set of active DOFs: 
%      - POINTSl: free (movable and optimizable) points.
%      - POINTSRpw: fixed (zero-weight) points.
%   2. While convergence is achieved and point count > target:
%      - Calls `MAKE1ZERO` to attempt removal of the least significant point.
%      - Updates weights and point positions to satisfy quadrature constraints.
%      - Records history (positions, weights, errors).
%
% REMARKS:
%   - `MAKE1ZERO` performs the Newton-Raphson correction and error monitoring.
%   - If a maximum number of points is specified in DATALOC.NumberOfCECMpoints,
%     the process stops when this number is reached.
%   - The output quadrature rule is guaranteed (within tolerance) to satisfy
%     the reduced basis moment conditions.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp1.mat')
end

CONVERGENCE = 1 ;
iter = 1;
VARC.POINTSl = find(wOLD ~=0);% These are the indexes of the points for which both weights and positions are unknowns
% (candidates for being eliminated)
VARC.POINTSRp = []; % These are the indexes of the points for which   positions   are constrained, but weights aren't
VARC.POINTSRpw = find(wOLD ==0); % These are the indexes of the points for which   weights and positions are constrained
DATALOC = DefaultField(DATALOC,'Include2ndStageIterations_PlotEvolutionWeights',1) ; % Plotting option
DATALOC = DefaultField(DATALOC,'VARIABLE_ITERATIONS_SECOND_STAGE',0) ; % Number of iterations is variable, see Eliminate1PointPROG.m

HISTORY.INDEXES_FIRST_STAGE = 1:length(HISTORY.POINTS_all) ; % For Plotting purposes
HISTORY.INDEXEX_SECOND_STAGE = {};% For Plotting purposes
HISTORY.CONTROL_POINTS = [] ; 
HISTORY.POINTS = HISTORY.POINTS_all ; 
HISTORY.WEIGHTS = HISTORY.WEIGHTS_all ; 

%HISTORY.CONTROL_POINTS =[] ;
NumberIterFirstStage = length(HISTORY.INDEXES_FIRST_STAGE) ;
CurrentITERglo = NumberIterFirstStage  ;
% ------------------------------------------------------
if isempty(DATALOC.NumberOfCECMpoints)
    DATALOC.NumberOfCECMpoints = 0  ; % 28-Apr-2024
end
while CONVERGENCE ==1  &&  length(find(wOLD>0)) >DATALOC.NumberOfCECMpoints
    DATALOC.iter = iter ;
    %  [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = ControlPointsAlgLARGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC) ;
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC ] = MAKE1ZERO(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO) ;
    if CONVERGENCE == 1
        
        [HISTORY,INDnneg,CurrentITERglo] = HistoryPointsUpdate2023(wNEW,xNEW,HISTORY,DATALOC,iter,CurrentITERglo) ;
        DATALOC.HISTORY_LOCAL.x = {} ;
        DATALOC.HISTORY_LOCAL.w  = {} ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
        
        if   length(find(wOLD>0))== DATALOC.NumberOfCECMpoints
            disp('Number of points equal to the number prescribed to the user. ..Exiting')
            break
        end
        
    end
end
%HISTORY.CONTROL_POINTS = DATALOC.HISTORY_LOCAL.ControlPoints ;

if iter == 1
    INDnneg = find(wOLD ~=0);
end

disp('------------------------------------------------------------------------------------------------')
disp(['SECOND STAGE: Integration rule with ',num2str(length(INDnneg)),' of ',num2str(length(wOLD)),' points'])
disp('------------------------------------------------------------------------------------------------')
 
