function [POLYINFO,ELEMENTS_CONTAINING_xNEW,PHIk_y] = LocationPointsInElement_noCN(xNEW,VAR_SMOOTH_FE,POLYINFO,PHIk_y)

npoints = size(xNEW,1) ; % Number of points at which to evaluate the function

[INDEXES_NEAR,IND_POLYG] = FindClosestNodes(xNEW,VAR_SMOOTH_FE) ; 
ELEMENTS_CONTAINING_xNEW = zeros((npoints),1) ; % Indices of the elements containing the points
inew = 1;
POLYINFO.ListPointsOutside = [] ; % Lists points outside
while inew <=npoints
    xLOC = xNEW(inew,:); % Coordinate of the point at which the function is to be evaluated
    INDnear = INDEXES_NEAR(inew) ;  % INDEX Nearest FE mesh  node to xLOC
    % Searching for the element containing xLOC --> elemCONTAINER
    [elemCONTAINER,POLYINFO ]= WhichElementInside2023(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG,POLYINFO,inew) ;  
        
    if isempty(elemCONTAINER)
       % Point outside the domain 
        PHIk_y = [];
        ELEMENTS_CONTAINING_xNEW(inew ) = 0 ;
        POLYINFO.ListPointsOutside = [POLYINFO.ListPointsOutside;inew] ; 
     
    else        
        ELEMENTS_CONTAINING_xNEW(inew ) = elemCONTAINER ;  % Element containing the point        
    end
    inew = inew +1 ;
end