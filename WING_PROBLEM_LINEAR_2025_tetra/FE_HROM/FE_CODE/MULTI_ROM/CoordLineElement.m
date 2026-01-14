function [xC,yC,L_total,CNb] = CoordLineElement(NODES_LINESall,MESH2D,iline)


 NODESlocALL = NODES_LINESall{iline} ; % Nodes (midside and vertices) formed by this line
    % Boundary elements corresponding to these nodes
    [CNb setBelem]= ElemBnd(MESH2D.CNbALL,NODESlocALL) ;  % Including midside nodes
    
    % Coordinates of all points 
    NODESall = unique(CNb(:)) ; 
    COORall = MESH2D.COORall(NODESall,:); 
    % Centroid 
    xC = 0.5*(max(COORall(:,1))+min(COORall(:,1)) ); 
    yC = 0.5*(max(COORall(:,2))+min(COORall(:,2)) ); 
    % Total Length 
    dX =  max(COORall(:,1))-min(COORall(:,1))  ; 
    dY =  max(COORall(:,2))-min(COORall(:,2))  ;
    L_total = norm([dX dY]) ;