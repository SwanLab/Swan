function HistoryPointsALL2D(AUXVAR,MESH,VAR_SMOOTH_FE)
%--------------------------------------------------------------------------
% function HistoryPointsALL2D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% PURPOSE:
%   Visualizes the evolution of cubature points during the CECM process
%   in 2D domains. This function routes plotting to appropriate graphical
%   routines (`HistoryPointsDECM2D`, `HistoryPointsPlot2D`) depending on
%   user-defined options such as video generation or iterative snapshot
%   visualization.
%
%   The visualization helps track the movement, selection, and eventual
%   elimination of integration points during the optimization phases of
%   empirical cubature schemes.
%
% INPUTS:
%   - AUXVAR          : Structure containing configuration flags, colors,
%                       marker styles, and history of point evolution.
%                       Fields used include:
%                         * DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY
%                         * DATA_from_MAIN.MAKE_VIDEO_POINTS
%                         * DATAOUTdecm.HistoryPoints
%                         * ShowLocationPointsBASE
%                         * ColorToBeEliminated
%                         * MarkerToBeEliminated
%                         * LineWidth
%                         * etc.
%
%   - MESH            : Mesh structure including nodal coordinates (MESH.COOR)
%                       and possibly connectivity (MESH.CN) for plotting.
%
%   - VAR_SMOOTH_FE   : Structure containing smoothed FE data or auxiliary
%                       variables used for graphical overlay or legend details.
%
% FUNCTIONAL FLOW:
%   - Ensures all necessary fields in `AUXVAR` are initialized using
%     `DefaultField` utility.
%   - If enabled, calls:
%       * `HistoryPointsDECM2D` for video or dynamic visualization of point evolution.
%       * `HistoryPointsPlot2D` for static frame visualization.
%   - If DECM visualization is not enabled, it defaults to `HistoryPointsPlot2D`.
%
% OUTPUT:
%   - This function produces figures and/or animations but has no direct output.
%
% DEPENDENCIES:
%   - HistoryPointsPlot2D
%   - HistoryPointsDECM2D
%
% EXAMPLE USAGE:
%   HistoryPointsALL2D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% NOTES:
%   - This function assumes a 2D spatial domain based on the mesh dimension.
%   - Intended to be called from within `HistoryPointsPlot` after automatic
%     domain dimensionality detection.
%
% SEE ALSO:
%   HistoryPointsPlot, HistoryPointsPlot2D, HistoryPointsDECM2D,
%   HistoryPointsALL1D, HistoryPointsALL3D
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
    close all
    AUXVAR.ShowLocationPointsBASE =1 ; %
    AUXVAR.AmplificationFactorDomain =0  ;
    AUXVAR.ASPECT_RATIO =[1,1,0.5]  ;
    AUXVAR.LineWidth =5 ;
    AUXVAR.ColorLines  = [0,0,0] ;
    AUXVAR.ColorToBeEliminated  = [1,0,0] ;
    AUXVAR.MarkerToBeEliminated  = 'o' ;
    AUXVAR.MarkerSizeToBeEliminated  = 6 ;
    AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = 1;
    AUXVAR.DATA_from_MAIN.DelayBetweenFrames =0.2;
    AUXVAR.DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY =1;
    
        AUXVAR.DATA_from_MAIN.HIDE_LEGEND =1;

    
end

AUXVAR = DefaultField(AUXVAR,'DATAOUTdecm',[]) ;
AUXVAR.DATAOUTdecm = DefaultField(AUXVAR.DATAOUTdecm,'HistoryPoints',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY',1) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'MAKE_VIDEO_POINTS',0) ;
AUXVAR = DefaultField(AUXVAR,'ShowLocationPointsBASE',1) ;
AUXVAR = DefaultField(AUXVAR,'AmplificationFactorDomain',0.2) ;
AUXVAR = DefaultField(AUXVAR,'LineWidth',4) ;
AUXVAR = DefaultField(AUXVAR,'ColorLines',[0,0,0]) ;
AUXVAR = DefaultField(AUXVAR,'ASPECT_RATIO',[1,1,0.5]) ;
AUXVAR = DefaultField(AUXVAR,'ColorToBeEliminated',[1,0,0]) ;
AUXVAR = DefaultField(AUXVAR,'MarkerToBeEliminated','x') ;
AUXVAR = DefaultField(AUXVAR,'MarkerSizeToBeEliminated',10) ;
 
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'HIDE_LEGEND',0) ;






if AUXVAR.DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY == 1 && ~isempty( AUXVAR.DATAOUTdecm.HistoryPoints)
    if AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS == 1
        HistoryPointsDECM2D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
        %  HistoryPointsPlot1D(AUXVAR,MESH) ;
    else
        HistoryPointsDECM2D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
        HistoryPointsPlot2D(AUXVAR,MESH) ;
    end
    
else
    HistoryPointsPlot2D(AUXVAR,MESH) ;
end
