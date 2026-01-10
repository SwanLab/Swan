function HistoryPointsALL3D(AUXVAR,MESH,VAR_SMOOTH_FE)
%--------------------------------------------------------------------------
% function HistoryPointsALL3D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% PURPOSE:
%   Visualizes the evolution of cubature points during the Continuous
%   Empirical Cubature Method (CECM) in 3D domains. Depending on the
%   configuration set in AUXVAR, the function provides either animated or
%   static plots illustrating how the cubature rule adapts iteratively.
%
%   The goal is to give visual feedback on the selection, movement, and
%   elimination of quadrature points throughout the cubature optimization
%   process, especially in high-dimensional domains.
%
% INPUTS:
%   - AUXVAR          : Structure with configuration options and plotting
%                       styles. Expected subfields:
%                         * DATAOUTdecm.HistoryPoints
%                         * DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY
%                         * DATA_from_MAIN.MAKE_VIDEO_POINTS
%
%   - MESH            : Structure containing the coordinates of the full mesh
%                       (MESH.COOR) and potentially connectivity (MESH.CN),
%                       used for drawing the domain and context geometry.
%
%   - VAR_SMOOTH_FE   : Structure containing auxiliary smoothed data used
%                       for enhanced visualization (e.g., basis functions or
%                       domain overlays).
%
% FUNCTIONAL FLOW:
%   1. Initializes optional fields in `AUXVAR` using `DefaultField`.
%   2. Based on `SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY`, chooses between:
%       a) `HistoryPointsDECM3D` for animation or incremental display.
%       b) `HistoryPointsPlot3D` for final snapshot visualization.
%   3. If no DECM history exists or visualization is disabled, only static
%      plotting is performed.
%
% OUTPUT:
%   - No direct output. The function generates visual plots or videos as side effects.
%
% DEPENDENCIES:
%   - HistoryPointsDECM3D
%   - HistoryPointsPlot3D
%   - DefaultField (utility function for field initialization)
%
% EXAMPLE USAGE:
%   HistoryPointsALL3D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% NOTES:
%   - This routine is automatically triggered from `HistoryPointsPlot` if the
%     mesh dimension is 3.
%   - Animation is enabled via `MAKE_VIDEO_POINTS`.
%   - `HistoryPointsDECM3D` and `HistoryPointsPlot3D` must be defined elsewhere
%     for this function to operate.
%
% SEE ALSO:
%   HistoryPointsPlot, HistoryPointsALL2D, HistoryPointsALL1D,
%   HistoryPointsPlot3D, HistoryPointsDECM3D
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp2.mat')
    AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = 0; 
end


AUXVAR = DefaultField(AUXVAR,'DATAOUTdecm',[]) ;
AUXVAR.DATAOUTdecm = DefaultField(AUXVAR.DATAOUTdecm,'HistoryPoints',[]) ;
AUXVAR = DefaultField(AUXVAR,'DATA_from_MAIN',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY',1) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'MAKE_VIDEO_POINTS',0) ;

if AUXVAR.DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY == 1 && ~isempty( AUXVAR.DATAOUTdecm.HistoryPoints)
    if AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS == 1
        HistoryPointsDECM3D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
    else
        HistoryPointsDECM3D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
        HistoryPointsPlot3D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
    end
else
    HistoryPointsPlot3D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
end
