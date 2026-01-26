function [HISTORY,INDnneg,CurrentITERglo] = HistoryPointsUpdate2023(wNEW,xNEW,HISTORY,DATALOC,iter,CurrentITERglo)
%--------------------------------------------------------------------------
% function [HISTORY, INDnneg, CurrentITERglo] = ...
%     HistoryPointsUpdate2023(wNEW, xNEW, HISTORY, DATALOC, iter, CurrentITERglo)
%
% PURPOSE:
%   Updates the history structure with the latest integration rule 
%   (points and weights) obtained during a second-stage iteration 
%   of the Continuous Empirical Cubature Method (CECM).
%
% INPUTS:
%   - wNEW             : [n x 1] Updated vector of integration weights.
%   - xNEW             : [n x ndim] Updated coordinates of integration points.
%   - HISTORY          : Structure storing all iterations of points/weights.
%   - DATALOC          : Control structure with plotting and local history flags.
%   - iter             : Current iteration number of the second stage.
%   - CurrentITERglo   : Global iteration counter across all stages.
%
% OUTPUTS:
%   - HISTORY          : Updated history structure, including full trace of:
%                        * POINTS{•}, WEIGHTS{•} (compressed history)
%                        * POINTS_all{•}, WEIGHTS_all{•} (full history)
%                        * INDEXEX_SECOND_STAGE{•} (for plotting sequences)
%                        * CONTROL_POINTS (if tracking evolution)
%                        * ISALLPOSITIVE (flags if all weights are positive)
%   - INDnneg          : Indices of points with non-zero weights (support of rule).
%   - CurrentITERglo   : Updated global iteration counter.
%
% NOTES:
%   - If `DATALOC.Include2ndStageIterations_PlotEvolutionWeights` is disabled,
%     only the current (xNEW, wNEW) is stored in full history.
%   - If enabled, intermediate Newton-Raphson iterations stored in 
%     `DATALOC.HISTORY_LOCAL.x` and `w` are also included.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------

INDnneg = find(wNEW ~=0);
 
 HISTORY.POINTS{end+1} = xNEW ;
    HISTORY.WEIGHTS{end+1} = wNEW  ;

if  DATALOC.Include2ndStageIterations_PlotEvolutionWeights == 0
    HISTORY.POINTS_all{end+1} = xNEW ;
    HISTORY.WEIGHTS_all{end+1} = wNEW  ;
else
    HISTORY.INDEXEX_SECOND_STAGE{iter} = (CurrentITERglo+1):(CurrentITERglo+size(DATALOC.HISTORY_LOCAL.x,2));
    CurrentITERglo = CurrentITERglo+size(DATALOC.HISTORY_LOCAL.x,2) ;
    for inewentries = 1:length(DATALOC.HISTORY_LOCAL.x)
        HISTORY.POINTS_all{end+1} =DATALOC.HISTORY_LOCAL.x{inewentries} ;
        HISTORY.WEIGHTS_all{end+1} = DATALOC.HISTORY_LOCAL.w{inewentries} ;
    end
    HISTORY.CONTROL_POINTS(end+1) = DATALOC.HISTORY_LOCAL.ControlPoints ;

  %   ; 

end


HISTORY.ISALLPOSITIVE(end+1) = all(wNEW(INDnneg)>0) ;
