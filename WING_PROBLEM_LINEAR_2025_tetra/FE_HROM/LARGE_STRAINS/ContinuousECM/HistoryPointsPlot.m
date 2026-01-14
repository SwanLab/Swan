function HistoryPointsPlot(AUXVAR,MESH,VAR_SMOOTH_FE)
%--------------------------------------------------------------------------
% function HistoryPointsPlot(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% PURPOSE:
%   Visualizes the evolution of selected integration (cubature) points
%   throughout the CECM or HROM procedure. It routes the visualization 
%   to specialized functions for 1D, 2D, or 3D domains depending on the
%   spatial dimension of the mesh.
%
%   The function is typically used after running the sparsification
%   and optimization routines, and provides insight into the history
%   of weight evolution, point movement, and selection during the 
%   cubature compression.
%
% INPUTS:
%   - AUXVAR         : Structure containing auxiliary variables and
%                      history data of the CECM iterations. Required fields:
%                      * DATALOC: contains plotting parameters and flags.
%                      * HISTORY: contains .POINTS_all and .WEIGHTS_all, 
%                                  which track the evolution of points/weights.
%
%   - MESH           : Structure with mesh information.
%                      Fields used:
%                        * COOR : [Nnodes x ndim] array of nodal coordinates.
%
%   - VAR_SMOOTH_FE  : Structure containing smoothed quantities from FE mesh
%                      (e.g., coordinates, connectivity, types) used for 
%                      visual enhancement and plotting.
%
% OUTPUTS:
%   No output arguments. This function generates plots directly by calling:
%     * HistoryPointsALL1D()  – for 1D problems
%     * HistoryPointsALL2D()  – for 2D problems
%     * HistoryPointsALL3D()  – for 3D problems
%
% NOTES:
%   - Assumes that the dimensionality of the problem is inferred from
%     the number of columns in `MESH.COOR`.
%   - Function behavior is modular and easily extended by adapting or
%     implementing the corresponding HistoryPointsALL*D routines.
%
% EXAMPLE USAGE:
%   HistoryPointsPlot(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% SEE ALSO:
%   HistoryPointsALL1D, HistoryPointsALL2D, HistoryPointsALL3D
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------

% Patterned after /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/10_HOMOG_3D_1T/WeightsPlotEvolution.m
% more specifically:
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/10_HOMOG_3D_1T/WeightPlot2D.m

if nargin == 0
    load('tmp1.mat')
    close all
    
end

if size(MESH.COOR,2) == 1
 %   HistoryPointsPlot1D(AUXVAR,MESH) ;
    HistoryPointsALL1D(AUXVAR,MESH,VAR_SMOOTH_FE)
elseif size(MESH.COOR,2) == 2
    HistoryPointsALL2D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
elseif size(MESH.COOR,2) == 3
    
    HistoryPointsALL3D(AUXVAR,MESH,VAR_SMOOTH_FE) ; 
    
end