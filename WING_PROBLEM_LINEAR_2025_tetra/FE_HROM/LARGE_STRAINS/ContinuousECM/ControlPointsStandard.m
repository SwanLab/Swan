function  [xGOOD,wGOOD,DATALOC,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPointsStandard(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT)

CONVERGENCE = 1 ; 
 iter = 1;
 wNEW = [] ;  
 INDnneg = 1:length(wOLD) ;
 xGOOD = xOLD ; wGOOD = wOLD ; 
% ------------------------------------------------------
while CONVERGENCE ==1  %&& length(xOLD)>1
    DATALOC.iter = iter ;     
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = ControlPointsAlgLARGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC) ;
    if CONVERGENCE == 1        
        INDnneg = find(wNEW ~=0);     
        ELEMENTS_TO_PLOT = unique([ELEMENTS_TO_PLOT;POLYINFO.ELEMENTS_CONTAINING_xNEW(INDnneg) ] ) ;
        HISTORY.POINTS{end+1} = xNEW(INDnneg,:) ;      
        HISTORY.WEIGHTS{end+1} = wNEW(INDnneg,:) ;
        HISTORY.POINTS_all{end+1} = xNEW ;      
        HISTORY.WEIGHTS_all{end+1} = wNEW  ;
       
        ALLPOSI = all(wNEW(INDnneg)>0)  ; 
        if ALLPOSI
            xGOOD = xNEW ; 
            wGOOD = wNEW ; 
        end
         HISTORY.ISALLPOSITIVE(end+1) = ALLPOSI ; 
        HISTORY.ELEMENTS_CONTAINING_POINTS{end+1} =  POLYINFO.ELEMENTS_CONTAINING_xNEW(INDnneg)   ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
        if length(INDnneg) == 1
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

