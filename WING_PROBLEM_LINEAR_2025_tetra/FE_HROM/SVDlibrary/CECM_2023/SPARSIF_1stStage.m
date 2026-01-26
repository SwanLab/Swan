function  [xGOOD,wGOOD,DATALOC,POLYINFO,VARC,HISTORY] = SPARSIF_1stStage(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY)
%--------------------------------------------------------------------------
% function [xGOOD,wGOOD,DATALOC,POLYINFO,VARC,HISTORY] = ...
%              SPARSIF_1stStage(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY)
%
% PURPOSE:
%   Implements the **first stage** of the Continuous Empirical Cubature Method (CECM),
%   based on iterative **point elimination** using Newton-type updates.
%   In each iteration, a single cubature point is forced to zero weight while 
%   attempting to preserve exact integration of the reduced snapshot basis.
%
%   This stage yields a sparser integration rule with a reduced number of
%   integration points while maintaining integration accuracy.
%
% INPUT:
%   - xOLD           : [n x ndim] Coordinates of the initial cubature points (typically DECM).
%   - wOLD           : [n x 1] Initial positive weights (from DECM).
%   - b              : [r x 1] Exact integrals of reduced basis snapshots: b = PHI' * W.
%   - DATALOC        : Struct with parameters for Newton-Raphson elimination strategy.
%                      Includes:
%                         .TOL_NewtonRaphson_EliminationPoints
%                         .MaxIterationsNR_ElimPoints
%                         .NumberOfCECMpoints (optional target number of points)
%   - VAR_SMOOTH_FE  : Struct containing geometry, basis, connectivity, etc.
%   - VARC           : Structure defining which point DOFs are active/inactive in optimization.
%   - HISTORY        : Struct for tracking all intermediate states during reduction.
%
% OUTPUT:
%   - xGOOD          : [m x ndim] Coordinates of final retained points with positive weights.
%   - wGOOD          : [m x 1] Final weights of retained integration points.
%   - DATALOC        : Updated data structure.
%   - POLYINFO       : Updated polynomial evaluation structure (e.g., element mapping).
%   - VARC           : Updated design variable structure.
%   - HISTORY        : Struct storing the evolution of points, weights, and positivity status.
%
% METHOD:
%   - The function calls `MAKE1ZERO_1STEP` iteratively to eliminate one point per step.
%   - It keeps only the points with nonzero weights, updating their weights to
%     preserve the moment-matching condition PHI' * w ≈ b.
%   - It stops when the number of points reaches the prescribed target, or when
%     only one point remains.
%
% TERMINATION CONDITIONS:
%   - The user-defined number of points `NumberOfCECMpoints` is reached.
%   - Only one point remains (minimal set).
%   - Maximum Newton iterations or residual tolerance is satisfied.
%
% NOTE:
%   - The output may still be refined in a **second stage** (SPARSIF) if enabled.
%   - If all remaining weights are positive, the state is marked as "all positive".
%
% REFERENCE:
%   J.A. Hernández et al., "Continuous Empirical Cubature Method for Data Compression
%   in Computational Mechanics", Computer Methods in Applied Mechanics and Engineering, 2024.
%
% AUTHOR:
%   Joaquín A. Hernández (UPC-CIMNE), April 2023.
%--------------------------------------------------------------------------

% First stage of reduction (direct removal of weights)
% .....................................................
if nargin == 0
    load('tmp.mat')
end

CONVERGENCE = 1 ;
iter = 1;
wNEW = [] ;
INDnneg = 1:length(wOLD) ;
xGOOD = xOLD ; wGOOD = wOLD ;
POLYINFO.setElements = VAR_SMOOTH_FE.setElements; %  Set of elements containing xNEW (initial guess)
POLYINFO.TriangulationDelaunay = cell(size(VAR_SMOOTH_FE.CN,1),1) ;   % Triang. For Delaunay 3D search 
DATALOC = DefaultField(DATALOC,'TOL_NewtonRaphson_EliminationPoints',1e-8) ;
DATALOC = DefaultField(DATALOC,'MaxIterationsNR_ElimPoints',40) ;

% ------------------------------------------------------
while CONVERGENCE ==1  %&& length(xOLD)>1
    DATALOC.iter = iter ;
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = MAKE1ZERO_1STEP(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO) ;
    if CONVERGENCE == 1
        INDnneg = find(wNEW ~=0);
        HISTORY.POINTS_all{end+1} = xNEW ;
        HISTORY.WEIGHTS_all{end+1} = wNEW  ;        
        ALLPOSI = all(wNEW(INDnneg)>0)  ;
        if ALLPOSI
            xGOOD = xNEW ;
            wGOOD = wNEW ;
        end
        HISTORY.ISALLPOSITIVE(end+1) = ALLPOSI ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
        if length(INDnneg) == 1 
            break
        end
        
        if length(INDnneg) ==   DATALOC.NumberOfCECMpoints
            disp(['Number of CECM points equal to the number prescribed by the user...'])
            break
        end
        
    end
end

if ~isempty(wNEW)
    disp('------------------------------------------------------------------------------------------------')
    disp(['FIRST STAGE: Integration rule with ',num2str(length(INDnneg)),' of ',num2str(length(wNEW)),' points'])
    disp('------------------------------------------------------------------------------------------------')
else
    error(['Convergence not achieved; try to decrease the tolerancd for the Newton-Raphson'])
end

