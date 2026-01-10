function [inELEM,onELEM] = IsInsideALL(xLOC,VAR_SMOOTH_FE,elemLOC,IND_POLYG) 


 if length(xLOC) ==1
            % 1D problem
            elemCN = VAR_SMOOTH_FE.CN(elemLOC,:) ;
            COORelemLOC = VAR_SMOOTH_FE.COOR(elemCN,:) ; 
             
            inELEM = 0 ; onELEM = 0 ; 
            if  xLOC >= COORelemLOC(1) && xLOC <= COORelemLOC(2)
                inELEM = 1; 
            end          
        else
        [inELEM,onELEM] = IsInsideGeneral(xLOC,VAR_SMOOTH_FE.COOR,VAR_SMOOTH_FE.CN,elemLOC,IND_POLYG) ;
        end