function [NODESfaces,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
    = RotationModes(idomain,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
    COORrve)
% See DomainDecom_SVD.m
% -----------------------
%dbstop('7')
if nargin == 0
    load('tmp.mat')
end
ndim = size(COORabs,2) ;
xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;   % Boundaries of domain i (square)
xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
if ndim == 2
    zmin = [] ;   zmax = [] ;
else
    zmin = min(COORabs(:,3)) ; zmin = zmin(1) ;
    zmax = max(COORabs(:,3)) ; zmax = zmax(1) ;
end
% Boundaries of each cell
% Determination of "REFERENCE POINT"
switch  DATA.TypeUnitCell
    case {'HEXAG_2D_SQUARE','HEXAGONAL'}
        [NODESfaces{idomain} NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
        % REFERENCE POINT  (TO MEASURE COORDINATES, AS WELL AS DISPLACEMENTS)
        if ndim==3
            [NODESln  PLANESlines]= DetermineLinesPeriodicHEXAbcs(NODESfaces{idomain}) ;
            nnode = size(COORabs,1) ;
            lbott = NODESln{3} ;
            [xminREF NODEREF] = min(COORabs(lbott,2)) ;
            NODEREF = lbott(NODEREF) ;
            
            % Node 2
            [xminREF NODEREF2] = max(COORabs(lbott,2)) ;
            NODEREF2 = lbott(NODEREF2) ;
            % Node 3
            lbott = NODESln{5} ;
            [xminREF NODEREF3] = max(COORabs(lbott,2)) ;
            NODEREF3 = lbott(NODEREF3) ;
            
            
            
        else
            
            % BOTTOM NODES, leftmost
            lbott = NODESfaces{idomain}{2} ;
            [xminREF NODEREF] = min(COORabs(lbott,1)) ;
            NODEREF = lbott(NODEREF) ;
            % Reference point 2  (top)
            lbott = NODESfaces{idomain}{4 } ;
            [xminREF NODEREF2] = min(COORabs(lbott,1)) ;
            NODEREF2 = lbott(NODEREF2) ;
        end
        
        
    otherwise
        error('Option not implemented')
end
% Displacements




if ndim == 3
    % Now coordinates and displacements are referred to the reference point
    
    
    % Reference node
    % ------------------
    
    
    RIGID_BODY_MOTION = repmat(dRVE{idomain}(NODEREF,:),nnode,1) ;
    CoordinatesChange = repmat(COORabs(NODEREF,:),nnode,1);  % Change of coordinates
    COORrve{idomain} = COORabs - CoordinatesChange ;
    dRVE{idomain} =  dRVE{idomain} -RIGID_BODY_MOTION  ;  % Avoid translations
    
    
    
    % dbstop('81')
    PointsRef = [NODEREF2';NODEREF3'] ;
    dREF = reshape(dRVE{idomain}(PointsRef,:)',[],1) ;
    %     xREF = reshape(COORrve{idomain}(PointsRef,:)',[],1) ;
    %     Omega = xREF\dREF ;
    % disp('Borrar')
    % d23 = dRVE{idomain}(PointsRef,:)' ;
    % x23 = COORrve{idomain}(PointsRef,:)' ;
    % OmegaEst = x23/d23  ;
    %
    % disp('End borrar')
    %
    
    Omega =[] ;
    for ipoint = 1:length(PointsRef)
        c1 = COORrve{idomain}(PointsRef(ipoint),:) ;
        OmegaLOC = [0 -c1(3) c1(2) ; c1(3) 0 -c1(1); -c1(2) c1(1) 0] ;
        Omega = [Omega; OmegaLOC] ;
    end
    %[~,idx]=licols(Omega') ;
    idx = [3 4 6] ;
    
    Omega = Omega(idx,:) ;
    
    
    ROTATIONS = Omega\dREF(idx) ;
    MODES_ROT = zeros(ndim*nnode,3) ;
    MODES_ROT(2:3:end,1) =  COORrve{idomain}(:,3) ;
    MODES_ROT(3:3:end,1) =  -COORrve{idomain}(:,2) ;
    
    MODES_ROT(1:3:end,2) =  -COORrve{idomain}(:,3) ;
    MODES_ROT(3:3:end,2) =  COORrve{idomain}(:,1) ;
    
    MODES_ROT(1:3:end,3) =  COORrve{idomain}(:,2) ;
    MODES_ROT(2:3:end,3) =  -COORrve{idomain}(:,1) ;
    RIGID_BODY_MOTION_rot =     ROTATIONS(1)*MODES_ROT(:,1) +  ROTATIONS(2)*MODES_ROT(:,2)  +ROTATIONS(3)*MODES_ROT(:,3) ;
    RIGID_BODY_MOTION_rot = reshape(RIGID_BODY_MOTION_rot,ndim,[])' ;
    
    dRVE{idomain} =  dRVE{idomain} -RIGID_BODY_MOTION_rot  ;
    RIGID_BODY_MOTION = RIGID_BODY_MOTION+RIGID_BODY_MOTION_rot ;
    
    % Matrix rotation ---skew matrix
    %     R3 = [0 -dX(3) dX(2); dX(3) 0 -dX(1); dX(2)  dX(1)  0] ;
    %     dispREF3 = dRVE{idomain}(NODEREF3,:) ;
    %     R =[R2; R3 ];
    %     b = [dispREF2'; dispREF3'];
    %     w = R\b ;
    %
    %     R2*w-dispREF2'
    
else
    
    
    
    %   Reference point 3  (top)
    % Now coordinates and displacements are referred to the reference point
    nnode = size(COORabs,1) ;
    RIGID_BODY_MOTION = repmat(dRVE{idomain}(NODEREF,:),nnode,1) ;
    CoordinatesChange = repmat(COORabs(NODEREF,:),nnode,1);  % Change of coordinates
    COORrve{idomain} = COORabs - CoordinatesChange ;
    dRVE{idomain} =  dRVE{idomain} -RIGID_BODY_MOTION  ;  % Avoid translations
    
    
    
    % Purging rotations
    dispHORIZONTAL2 = dRVE{idomain}(NODEREF2,1) ;  % Horizontal displacement 2nd reference point
    distanceBETWEEN_refP = (COORrve{idomain}(NODEREF2,2) - COORrve{idomain}(NODEREF,2))  ;
    ANGLE_ROTATION = dispHORIZONTAL2/distanceBETWEEN_refP;
    
    
    
    disp_HORIZONTAL = ANGLE_ROTATION*COORrve{idomain}(:,2) ;
    disp_VERTICAL = ANGLE_ROTATION*COORrve{idomain}(:,1) ;
    dRVE{end}(:,1) =  dRVE{idomain}(:,1) - disp_HORIZONTAL ;  % Avoid rotations
    dRVE{end}(:,2) =  dRVE{idomain}(:,2) + disp_VERTICAL ;  % Avoid rotations
    RIGID_BODY_MOTION(:,1) =    RIGID_BODY_MOTION(:,1) + disp_HORIZONTAL ;
    RIGID_BODY_MOTION(:,2) =    RIGID_BODY_MOTION(:,2) - disp_VERTICAL ;
end

%   RIGID_BODY_MOTION = reshape(RIGID_BODY_MOTION',[],1) ;

RIGID_BODY_MOTIONglo{idomain} = RIGID_BODY_MOTION ;
CoordinatesChangeGLO{idomain} = CoordinatesChange ;


