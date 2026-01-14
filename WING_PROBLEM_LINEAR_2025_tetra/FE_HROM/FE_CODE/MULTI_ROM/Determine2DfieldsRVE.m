function DISP2D= Determine2DfieldsRVE(a,MESH2D,ndim,MESH3D,DOFsKEEP,ndimGLO,DATAROM_glo,DATA_REFMESH)
if nargin == 0
    load('tmp2.mat')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a --> Displacement midside points of each quadrilateral element
% MESH2D --> Created in function Geometry2Dstructure.m
% Quadrilateral element determined by four planes x=xmin, x=xmax, y=ymin, y=ymax
% Midside points corresponds to columns 5 to 8 in CN, are numbered in the following order
% NODE 1 (5)  = plane x=xmin;  NODE 2 (6)= plane y=ymin, NODE 3 (7) = plane x=xmax, NODE 4 (8) = plane y=ymax
%
% % Corner points, on the other hand, are labelled as
% NODE 1 = planes (xmin,ymax), NODE 2 = planes (xmin,yin), NODE 3 = planes
% (xmax,ymin), NODE 4 = planes (xmax,ymax )
%
% The vector of displacement "a" corresponds to the generalized
% displacements of the centroid of the faces corresponding to its  midside points
% The location of the midside points may not coincide with the centroids of
% the faces
% This function computes, by interpolation, an estimation of the
% displacement of the remaining nodes (corner points)
% JAHO, 26-July-2018
% ---------------------------------------------------------------

%ndim = 3 ;
if ndim == 3
    nrigidbody = 6 ;
else
    nrigidbody = 3 ;
end

if ~isempty(DOFsKEEP)
    ndimMAX = max(ndimGLO) ;
    ndof = size(MESH2D.COOR,1)*ndimMAX ;
    aALL = zeros(ndof,1) ;
    aALL(DOFsKEEP) = a ;
    aALL = reshape(aALL,ndimMAX,[]) ;
else
    ndimGLO = ndimGLO(1);
    aALL = reshape(a,ndimGLO,[]) ;
end
%
if ~isstruct(DATAROM_glo{1}.BasisInt)
    % ************
    % OLD METHOD
    % ************
    a = aALL(1:nrigidbody,:) ;
    
    % ROTATIONS
    %----------
    %  DATA_REFMESH{1} = DefaultField(DATA_REFMESH{1},'RotationMatrixFace',[]) ;
    ROT_RELATIVE = [] ; ROT_RELATIVE_3 = [] ;
    if isfield(DATA_REFMESH{1},'RotationMatrixFace')
        faceROTATED = 3 ;
        ROT_RELATIVE = DATA_REFMESH{1}.RotationMatrixFace ;
        
        ROT_RELATIVE_3 = DATA_REFMESH{1}.RotationMatrixFace{faceROTATED} ;
    end
    
    
    %  [AAA,BBB ] = cellfun(@size,ROT_RELATIVE);
    
    if   ~isempty(ROT_RELATIVE_3)
        
        % We have to rotate each pINTF
        % First we have to determine the domain to which each node pertains
        % (This should be VECTORIZED in future versions)
        aNEW = zeros(size(a)) ;
        
        
        for inode = 1:size(aNEW,2)
            [eDOM,iNODEloc] =  find(inode ==  MESH2D.CN) ;
            
            eDOM = eDOM(1) ;
            iNODEloc = iNODEloc(1) ;
            Rot_GLO_DOMloc = MESH2D.rotDOM{eDOM} ;
            Rot_DOMloc_INTFloc = ROT_RELATIVE{iNODEloc}  ;
            if ~isempty(Rot_DOMloc_INTFloc)
                Rot_GLO_INTFloc = Rot_GLO_DOMloc*Rot_DOMloc_INTFloc ;
            else
                Rot_GLO_INTFloc = Rot_GLO_DOMloc ;
            end
            
            aNEW(1:3,inode) = Rot_GLO_INTFloc*a(1:3,inode) ;
            aNEW(4:6,inode) = Rot_GLO_INTFloc*a(4:6,inode) ;
            
        end
        a = aNEW ;
        
        
    end
    
    
else
    % not useful....
    % **************************************************
    % NEW METHOD, KINEMATICALLY CONSTRAINED , APR-2019
    % **************************************************
    % a --> contains the information of the MASTER dofs
    % This should be vectorized
    a = zeros(nrigidbody,size(aALL,2)) ;
    Acomp = DATAROM_glo{1}.BasisInt.Acomp  ;
    IndicesRB = DATAROM_glo{1}.BasisInt.IndicesRB ;
    nDOFsFACE = DATAROM_glo{1}.BasisInt.nDOFsFACE ;
    DOFsP =  DATAROM_glo{1}.BasisInt.DOFsP ;
    DOFmP =  DATAROM_glo{1}.BasisInt.DOFmP ;
    for ielem = 1:size(MESH2D.CN,1)
        CNloc = MESH2D.CN(ielem,:) ;
        pMASTERe = aALL(:,CNloc) ;
        pSLAVEe = zeros(size(Acomp,1),1) ;
        pMASTEReCOL = [];
        iini = 1;
        for iface= 1:length(DATAROM_glo{1}.BasisInt.BasisINTFall_cell)
            nDOFloc = nDOFsFACE(iface) ;
            ifin = iini+nDOFloc-1 ;
            ENTRIES = iini:ifin ;
            pSLAVEe =  pSLAVEe  + Acomp(:,ENTRIES)*pMASTERe(1:nDOFloc,iface) ;
            pMASTEReCOL = [pMASTEReCOL;pMASTERe(1:nDOFloc,iface)] ;
            iini = ifin+1;
        end
        
        pINTFall = zeros(length(pSLAVEe)+length(pMASTERe),1) ;
        
        pINTFall(DOFmP) = pMASTEReCOL ;
        pINTFall(DOFsP) = pSLAVEe ;
        
        pINTFrb = pINTFall(IndicesRB) ;
        pINTFrb = reshape(pINTFrb,nrigidbody,[]);
        
        a(:,CNloc) = pINTFrb ;
        
    end
    
    
end






nnode = size(MESH2D.COORall,1) ; % Number of nodes (quadratic mesh)
DISP2D = zeros(nnode,ndim) ; % Required output (displacement all nodes quadratic mesh )

NODES_MIDSIDE = MESH2D.NODESmid ;  % Global numbering of midside points
% Displacement centroids
% ----------------------
% The displacements of the centroids are given by the vector "a"
% Now we wish to infer the displacements of the midpoints (t)
MIDPOINTS = MESH3D.RVES.DATA3D.COOR_MIDSIDE_FACE ;
CENTROIDS  = MESH3D.RVES.DATA3D.CENTRf ;
CORNERS  =  MESH3D.RVES.DATA3D.COOR_CORNER ;
for  iii = 1:length(MIDPOINTS)
    COORrel{iii} = MIDPOINTS{iii}(:)-CENTROIDS{iii}(:) ;
    R = ConstructBasisRigidBody(COORrel{iii}') ;
    NODESquadratic = MESH2D.CNall(:,iii+4) ;  % Face-iii, list of nodes QUADRATIC element xº
    NODES_loc = MESH2D.CN(:,iii) ;
    aMIDSIDE = a(:,NODES_loc) ;
    dLOC = R*aMIDSIDE ;
    DISP2D(NODESquadratic,:) = dLOC' ;
end
%
% ------------------------------
% DISPLACEMENTS CORNER POINTS (computed element-wise)
% ----------------------------
% --------------
imidside = {[4, 1],[1,2],[2,3],[3,4] };
d = cell(size(imidside)); % Displacements
NODES_VERTEX = cell(size(imidside)); % Displacements
for ipoint = 1:length(imidside) % Loop over corner points
    [d{ipoint},NODES_VERTEX{ipoint}]= EstimationDisplacementCornerGLO(MESH2D,ipoint,imidside{ipoint},a,...
        MIDPOINTS,CENTROIDS,CORNERS) ;
end

% Next we perform a loop over all elements
NODESvert = MESH2D.NODESvert ;
DISP_CORNERS = zeros(nnode,ndim) ;
NODES_VISITED = zeros(nnode,1) ;
for inode = 1:length(NODESvert)
    node = NODESvert(inode) ;
    for ivertices = 1:length(NODES_VERTEX)
        III =  find(NODES_VERTEX{ivertices}==node) ;
        if ~isempty(III)
            DISP_CORNERS(node,:) = DISP_CORNERS(node,:) + d{ivertices}(:,III)' ;
            NODES_VISITED(node) = NODES_VISITED(node) + +1;
        end
    end
end

for idim  =1:ndim
    DISP_CORNERS(NODESvert,idim) = DISP_CORNERS(NODESvert,idim)./ NODES_VISITED(NODESvert) ;
end


DISP2D(NODESvert,:) = DISP_CORNERS(NODESvert,:) ;

DISP2D = DISP2D' ;


end



function d = EstimationDisplacementCorner(MESH2D,a,imidside1,COORcorner)

% Generalized displacement midside points plane
NODES_mid = MESH2D.CN(:,imidside1) ;
aMIDSIDE = a(:,NODES_mid) ;
% Coordinates midside points
COORmids = MESH2D.COOR(NODES_mid(1),1) ;
% Relative coordinates
COORrel = COORmids-COORcorner ;
% Rigid body modes
R = ConstructBasisRigidBody(COORrel) ;
% Displacement estimated by superposing translation  and rotation
d = R*aMIDSIDE ;

end

function  [d,NODES]= EstimationDisplacementCornerGLO(MESH2D,ipoint,imidside,a,MIDPOINTS,CENTROIDS,CORNERS)

NODES = MESH2D.CNall(:,ipoint) ;  % Global numbering  corner points "ipoint"

% 1) ESTIMATION DISPLACEMENT USING INFORMATION  iface = imidside(1) ;
iface = imidside(1) ;
COORrel = CORNERS{ipoint}(:) -CENTROIDS{iface}(:) ;
NODES_mid = MESH2D.CN(:,iface) ;
aMIDSIDE = a(:,NODES_mid) ;
R = ConstructBasisRigidBody(COORrel') ;
% Displacement estimated by superposing translation  and rotation
d1 = R*aMIDSIDE ;
% 22) ESTIMATION DISPLACEMENT USING INFORMATION PLANE ymax
iface = imidside(2) ;
COORrel = CORNERS{ipoint}(:) -CENTROIDS{iface}(:) ;
NODES_mid = MESH2D.CN(:,iface) ;
aMIDSIDE = a(:,NODES_mid) ;
R = ConstructBasisRigidBody(COORrel') ;
% Displacement estimated by superposing translation  and rotation
d4 = R*aMIDSIDE ;
% Average
d = (d1+d4)/2 ;

end
