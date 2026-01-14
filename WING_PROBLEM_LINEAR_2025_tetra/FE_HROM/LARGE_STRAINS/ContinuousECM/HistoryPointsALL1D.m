function HistoryPointsALL1D(AUXVAR,MESH,VAR_SMOOTH_FE)


AUXVAR = DefaultField(AUXVAR,'DATAOUTdecm',[]) ;
AUXVAR.DATAOUTdecm = DefaultField(AUXVAR.DATAOUTdecm,'HistoryPoints',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY',1) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'MAKE_VIDEO_POINTS',0) ;
if AUXVAR.DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY == 1 && ~isempty( AUXVAR.DATAOUTdecm.HistoryPoints)
    if AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS == 1
        HistoryPointsDECM1D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
        %  HistoryPointsPlot1D(AUXVAR,MESH) ;
    else
        HistoryPointsDECM1D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
        HistoryPointsPlot1D(AUXVAR,MESH) ;
    end
    
else
    HistoryPointsPlot1D(AUXVAR,MESH) ;
end
