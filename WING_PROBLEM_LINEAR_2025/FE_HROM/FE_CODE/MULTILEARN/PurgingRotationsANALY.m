function [NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
    = PurgingRotationsANALY(idomain,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
    COORrve)
% See DomainDecom_SVD.m
% -----------------------
if nargin == 0
    load('tmp.mat')
end

xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;   % Boundaries of domain i (square)
xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
zmin = [] ;   zmax = [] ;

% Boundaries of each cell
% Determination of "REFERENCE POINT"
switch  DATA.TypeUnitCell
    case 'HEXAG_2D_SQUARE'
        [NODESfaces{idomain} NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
        % REFERENCE POINT  (TO MEASURE COORDINATES, AS WELL AS DISPLACEMENTS)
        % BOTTOM NODES, leftmost
        lbott = NODESfaces{idomain}{2} ;
        [xminREF NODEREF] = min(COORabs(lbott,1)) ;
        NODEREF = lbott(NODEREF) ;
        % Reference point 2  (top)
        lbott = NODESfaces{idomain}{4 } ;
        [xminREF NODEREF2] = min(COORabs(lbott,1)) ;
        NODEREF2 = lbott(NODEREF2) ;
        
    otherwise
        error('Option not implemented')
end
% --------------------------
% Rigid body modes  (2D)
% --------------------------
% Rotation around NODEREF
% Change coordinates
nnodeRVE = size(COORabs,1) ;  ndim = size(COORabs,2) ;
CoordinatesChange = repmat(COORabs(NODEREF,:),nnodeRVE,1);  %
COORrve{idomain} = COORabs - CoordinatesChange ;  % New coordinates

MODErb = zeros(ndim*nnodeRVE,3) ;
MODErb(1:ndim:end,1) = 1/sqrt(nnodeRVE);
MODErb(2:ndim:end,2) = 1/sqrt(nnodeRVE);
USE_centroid = 0 ;
if USE_centroid == 0
    
    RIGID_BODY_MOTION_fixed = repmat(dRVE{idomain}(NODEREF,:),nnodeRVE,1) ;
    dRVE{idomain} =  dRVE{idomain} -RIGID_BODY_MOTION_fixed  ;
    
    d = reshape(dRVE{idomain}',nnodeRVE*ndim,1) ;
    
    
    xxx = COORrve{idomain}(:,1); %-COORgeo(:,1) ;
    yyy = COORrve{idomain}(:,2); %-COORgeo(:,2) ;
    MODErb(1:ndim:end,3) = -yyy;
    MODErb(2:ndim:end,3) = xxx;
    MODErb(:,3) = MODErb(:,3)/norm(MODErb(:,3));
    % -------------------------------------------------------
    % Rigid body motion
    RIGID_BODY_MOTION = MODErb*(MODErb\d)  ;
    d =  d -RIGID_BODY_MOTION  ;  %
    
    dRVE{idomain} = reshape(d,ndim,nnodeRVE)' ;
    RIGID_BODY_MOTION = reshape(RIGID_BODY_MOTION,ndim,nnodeRVE)' ;
    RIGID_BODY_MOTION = RIGID_BODY_MOTION + RIGID_BODY_MOTION_fixed ;
    
else
    d = reshape(dRVE{idomain}',nnodeRVE*ndim,1) ;
    
    xGEO = sum(COORabs(:,1))/(nnodeRVE) ;
    yGEO = sum(COORabs(:,2))/(nnodeRVE) ;
    COORgeo = repmat([xGEO,yGEO],nnodeRVE,1) ;
    xxx = COORabs(:,1)-COORgeo(:,1) ;
    yyy = COORabs(:,2)-COORgeo(:,2) ;
    MODErb(1:ndim:end,3) = -yyy;
    MODErb(2:ndim:end,3) = xxx;
    MODErb(:,3) = MODErb(:,3)/norm( MODErb(:,3)) ;
    RIGID_BODY_MOTION = MODErb*(MODErb'*d)  ;
    d =  d -RIGID_BODY_MOTION  ;  %
    dRVE{idomain} = reshape(d,ndim,nnodeRVE)' ;
    RIGID_BODY_MOTION = reshape(RIGID_BODY_MOTION,ndim,nnodeRVE)' ;
    
end
% Deformation






RIGID_BODY_MOTIONglo{idomain} = RIGID_BODY_MOTION ;
CoordinatesChangeGLO{idomain} = CoordinatesChange ;
