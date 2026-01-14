function [ELEMENTS_TO_PLOT,HISTORY,INDnneg,CurrentITERglo] = HistoryPointsUpdate(wNEW,xNEW,ELEMENTS_TO_PLOT,HISTORY,DATALOC,POLYINFO,iter,CurrentITERglo)

INDnneg = find(wNEW ~=0);
ELEMENTS_TO_PLOT = unique([ELEMENTS_TO_PLOT;POLYINFO.ELEMENTS_CONTAINING_xNEW(INDnneg) ] ) ;
HISTORY.POINTS{end+1} = xNEW(INDnneg,:) ;
HISTORY.WEIGHTS{end+1} = wNEW(INDnneg,:) ;

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
HISTORY.ELEMENTS_CONTAINING_POINTS{end+1} =  POLYINFO.ELEMENTS_CONTAINING_xNEW(INDnneg)   ;