function VariousPLOTS_CECM(DATALOC,DATAecm,MESH,AUXVAR,VAR_SMOOTH_FE)
% -----------------------------------------
if nargin == 0
    load('tmp1.mat')
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
        if DATALOC.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS ==3
            PrintingGid_ECMpoints_LARGE(MESH,DATALOC,AUXVAR.HISTORY) ;
            %     else
            %         error('Option not implemented')
        end
        
    end
    
    
    
    if size(coorECM,2) == 2
        PlotPathPointsGauss2D(AUXVAR.ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,AUXVAR.PATHiniPOINTS,coorECM,xNEW) ;
    elseif size(coorECM,2) == 3
        PlotPathPointsGauss3D(AUXVAR.ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,AUXVAR.PATHiniPOINTS,coorECM,xNEW) ;
        %elseif size(coorECM,2) == 1
        %   PlotPathPointsGauss1D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xNEW) ;
    end
    
    
    
    
    
    
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



