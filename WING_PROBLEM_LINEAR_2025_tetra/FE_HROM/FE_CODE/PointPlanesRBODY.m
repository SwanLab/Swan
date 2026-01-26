function [NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COORabs,CNref,DATA)


%dbstop('5')
if nargin ==0
    load('tmp2.mat')
end
f1NOD = [] ;
f2NOD = [] ;
ndim = size(COORabs,2)  ;


xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;   % Boundaries of domain i (square)
xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
if ndim ==2
    zmin = [] ;   zmax = [] ;
    DATA = DefaultField(DATA,'TypeUnitCell','HEXAG_2D_SQUARE') ;
    
else
    zmin = min(COORabs(:,3)) ; zmin = zmin(1) ;
    zmax = max(COORabs(:,3)) ; zmax = zmax(1) ;
    DATA = DefaultField(DATA,'TypeUnitCell','HEXAHEDRA') ;
    
end


TOL = ChooseTolerance(CNref,COORabs) ;
DATA = DefaultField(DATA,'CalculateMasterSlaves',1) ;


switch  DATA.TypeUnitCell
    case 'HEXAG_2D_SQUARE'
        [NODESfaces, NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
        % REFERENCE POINT  (TO MEASURE COORDINATES, AS WELL AS DISPLACEMENTS)
        % BOTTOM NODES, leftmost
        lbott = NODESfaces{2} ;
        [xminREF NODEREF] = min(COORabs(lbott,1)) ;
        NODEREF = lbott(NODEREF) ;
        % Reference point 2  (top)
        lbott = NODESfaces{4 } ;
        [xminREF NODEREF2] = min(COORabs(lbott,1)) ;
        NODEREF2 = lbott(NODEREF2) ;
        %%%
        if  DATA.CalculateMasterSlaves ==1
            [f1NOD f2NOD] =MasterSlavesSets2Dbeam(NODESfaces,COORabs,TOL,NormalPlanes,DATA) ;
        end
    case 'HEXAHEDRA'
        [NODESfaces, NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
       
        
        [NODESln  PLANESlines]= DetermineLinesPeriodicHEXAbcs(NODESfaces) ;
        nnode = size(COORabs,1) ;
        lbott = NODESln{3} ;
        [xminREF NODEREF] = min(COORabs(lbott,2)) ;
        NODEREF = lbott(NODEREF) ;
        %
        %             % Node 2
        %             [xminREF NODEREF2] = max(COORabs(lbott,2)) ;
        %             NODEREF2 = lbott(NODEREF2) ;
        %             % Node 3
        %             lbott = NODESln{5} ;
        %             [xminREF NODEREF3] = max(COORabs(lbott,2)) ;
        %             NODEREF3 = lbott(NODEREF3) ;
        
        
        if DATA.CalculateMasterSlaves ==1
            [f1NOD f2NOD NODESfaces] =MasterSlavesSets3Dbeam(NODESfaces,COORabs,TOL,NormalPlanes,DATA) ;
        end
        
    otherwise
        error('Option not implemented')
end