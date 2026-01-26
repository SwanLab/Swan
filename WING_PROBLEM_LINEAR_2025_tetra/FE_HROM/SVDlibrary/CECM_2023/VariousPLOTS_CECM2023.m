function VariousPLOTS_CECM2023(DATALOC,DATAecm,MESH,AUXVAR,VAR_SMOOTH_FE)
%--------------------------------------------------------------------------
% function VariousPLOTS_CECM2023(DATALOC,DATAecm,MESH,AUXVAR,VAR_SMOOTH_FE)
%
% PURPOSE:
%   Generates a set of diagnostic and illustrative plots related to the
%   Continuous Empirical Cubature Method (CECM), including visualization
%   of basis functions, integration point weights before and after optimization,
%   and the evolution of selected points throughout the CECM procedure.
%
% INPUTS:
%   - DATALOC        : Structure with plotting options and global parameters.
%                      Fields used: SHOW_FIGURES, PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS.
%
%   - DATAecm        : Structure containing CECM results.
%                      Fields used: 
%                         * xCECM     : Final optimized cubature points
%                         * wCECM     : Corresponding weights
%                         * xDECM     : Initial ECM points before optimization
%                         * wDECM     : Initial weights
%                         * setElements : Elements selected during reduction
%
%   - MESH           : Reference mesh structure used for visualization.
%
%   - AUXVAR         : Structure for auxiliary variables.
%                      Fields used:
%                         * DATA_from_MAIN.BasisIntegrandToPlot_INDEXES : indexes of integrand basis functions to visualize
%                         * ELEMENTS_TO_PLOT, PATHiniPOINTS (optional for 2D/3D path visualization)
%
%   - VAR_SMOOTH_FE  : Structure containing smoothed FE quantities used
%                      for plotting and diagnostics.
%                      Fields used:
%                         * BasisIntegrand : [Npoints x Nmodes] basis matrix
%                         * COORg          : spatial coordinate vector
%
% OUTPUTS:
%   No output arguments. The function generates visual plots directly.
%
% NOTES:
%   - If 1D basis functions are present and `BasisIntegrandToPlot_INDEXES`
%     is defined, a plot of selected basis functions is shown.
%   - Bar plots compare initial and final weights for visual diagnosis.
%   - Calls external routines such as `HistoryPointsPlot` for visualizing
%     the evolution of selected cubature points.
%
% EXAMPLE USAGE:
%   VariousPLOTS_CECM2023(DATALOC,DATAecm,MESH,AUXVAR,VAR_SMOOTH_FE)
%
% SEE ALSO:
%   HistoryPointsPlot.m, PlotPathPointsGauss2D.m, 
%   EvaluateBasisFunctionDIRECTFIT2023.m
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024.
%--------------------------------------------------------------------------

% -----------------------------------------
if nargin == 0
    load('tmp.mat')
end

ELEMENTS_xNEW = DATAecm.setElements;

coorECM = DATAecm.xDECM ;
wINITIAL = DATAecm.wDECM ;
xNEW = DATAecm.xCECM ;
wNEW = DATAecm.wCECM ;



AUXVAR = DefaultField(AUXVAR,'DATA_from_MAIN',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'BasisIntegrandToPlot_INDEXES',[]) ;
if ~isempty(AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES) && size(xNEW,2) == 1
    figure(564)
    hold on
    xlabel('x')
    ylabel('Integrand basis')
    title(['Basis functions for the integrand (total number = ',num2str(size(VAR_SMOOTH_FE.BasisIntegrand,2)),')'])
    h = zeros(size(AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES))
    LLL = cell(size(h)) ;
    for iplot = 1:length(AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES)
        iplotLOC = AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES(iplot) ;
        h(iplot) =   plot(VAR_SMOOTH_FE.COORg,VAR_SMOOTH_FE.BasisIntegrand(:,iplotLOC)) ;
        LLL{iplot} =['mode = ',num2str(iplotLOC)]; 
    end
    legend(h,LLL)
end




DATALOC.SHOW_FIGURES = 1;
if DATALOC.SHOW_FIGURES==1
    
    
    
    if ~isempty(ELEMENTS_xNEW)
        disp('MESH CONTAINING THE NEW REDUCED SET OF POINTS');
        %   HYPERREDUCED_VARIABLES.WdomRED = wNEW ;
        %  HYPERREDUCED_VARIABLES.setElements = ELEMENTS_xNEW ;
        %  DATAIN.LABEL_NAME_PROJECT = DATA_GENGAUSS.NAMEFILE ;
        %     DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_REDUCED_GAUSS_IN_GID',0) ;
        %     DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'SCALE_FACTOR_FOR_PRINTING_ECM_POINTS',0.1) ;
        %
        %     if DATA_GENGAUSS.PLOT_REDUCED_GAUSS_IN_GID == 1
        %         DATALOC.PRINT_MESH_ECM_POINTS_AS_POINTS = 1;
        %         DATALOC.COORDINATES_POINTS_TO_PRINT = xNEW ;
        %         DATALOC.COORDINATES_POINTS_TO_PRINT = xNEW ;
        %         DATALOC.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS = DATA_GENGAUSS.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS ;
        %     end
        %
        %
        DATALOC = DefaultField(DATALOC,'PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS',3) ; % =    1; % 0--> NONE % 1 -->Like time ... % 2 = Separated meshes
        
        %  if DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 0
        %     DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HROMVAR) ;
        % elseif DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 1
        %    DATAIN = PrintingGid_ECMhistory(DATAIN,DATA_REFMESH,HROMVAR,HISTORY_ITERATIONS) ;
%         if DATALOC.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS ==3
%             PrintingGid_ECMpoints_LARGE(MESH,DATALOC,AUXVAR.HISTORY) ;
%             %     else
%             %         error('Option not implemented')
%         end
        
    end
    
    
%     
%     if size(coorECM,2) == 2
%         PlotPathPointsGauss2D(AUXVAR.ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,AUXVAR.PATHiniPOINTS,coorECM,xNEW) ;
%     elseif size(coorECM,2) == 3
%         PlotPathPointsGauss3D(AUXVAR.ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,AUXVAR.PATHiniPOINTS,coorECM,xNEW) ;
%         %elseif size(coorECM,2) == 1
%         %   PlotPathPointsGauss1D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xNEW) ;
%     end
%     
%     
    
    
    
    
end


if size(coorECM,2) >1
    figure(16)
    hold on
    bar(sort(wNEW))
    ylabel('Weights (after optimization)')
    figure(17)
    hold on
    bar(sort(wINITIAL))
    ylabel('Weights (before optimization)')
    
else
    
end
%diary off

disp('------------------------------------------------')
disp('PLotting history points (and weights)')
disp('-------------------------------------------')
AUXVAR.DATALOC  = DATALOC ;
HistoryPointsPlot(AUXVAR,MESH,VAR_SMOOTH_FE) ; 



