function [INDEXES_NEAR,IND_POLYG] = FindClosestNodes(xNEW,VAR_SMOOTH_FE)

if size(xNEW,2) == 1
    INDEXES_NEAR = zeros(size(xNEW)) ;
    
    for ipointsX = 1:length(xNEW)
        [dummy,INDEXES_NEAR(ipointsX)] =  min(abs(VAR_SMOOTH_FE.COOR - xNEW(ipointsX))) ;
    end
    IND_POLYG = [1,2] ;
    
else
    INDEXES_NEAR = nearestNeighbor(VAR_SMOOTH_FE.DELTRIANG, xNEW);
    IND_POLYG = VAR_SMOOTH_FE.IND_POLYG_ELEMENT  ;  % Local numbering of corner nodes

end