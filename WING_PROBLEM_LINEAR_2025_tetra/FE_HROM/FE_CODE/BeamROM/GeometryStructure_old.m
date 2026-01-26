function GeometryStructure_old(NAMEMESH_STRUCTURE,NAMEMESH_SLICE,NAME_project,DATAIN)
% See BeamROM.tex, IMPLEM1.tex
%
% INPUTS : NAMEMESH_STRUCTURE
%          NAMEMESH_SLICE
%          NAME_project
%
%
% Reading coordinate structure (1D, linear)
[coorINT,CNint,TypeElement1D,CNbINT,TypeElementB,MaterialType1D]=...
    ReadMeshFile(NAMEMESH_STRUCTURE,'READ_MATERIAL_COLUMN',1)  ;

% INITIALIZATIONS 
% -------------------
% Type of subdomain (slice or joint)
TYPE_SUBDOMAIN_default = cell(size(NAMEMESH_SLICE))  ; 
TYPE_SUBDOMAIN_default(:) ={'slice'} ;
DATAIN = DefaultField(DATAIN,'TYPE_SUBDOMAIN',TYPE_SUBDOMAIN_default) ;
% 
ntypes = length(unique(MaterialType1D)) ; % number of distinct subdomains. 
if ntypes ~= length(NAMEMESH_SLICE); error('Check types of 1D element ') ; end
SLICEglo = cell(ntypes,1) ; % Cell of structures that will contain the geometric properties of the slices/joints
rotMATglo = cell(ntypes,1) ;% Cell of structures for storing the rotation matrices
% ---------------------

% Loop over number of slice/joint types  
 for itype = 1:ntypes
    % -----------------------------------------------
    % Reading mesh and geometric entities slice NAMEMESH_SLICE
    % Output 'NAME_SLICE', 'COOR','CN','TypeElement','CNb','TypeElementB','MaterialType'
    % 'NODES_face1','NODES_face2','CENTRf1','CENTRf2' (for slices)
    switch DATAIN.TYPE_SUBDOMAIN{itype}
        case   'slice'
            SLICE =  GeometrySlice(NAMEMESH_SLICE{itype},DATAIN) ;
        case 'joint'
             CheckMeshJoint(NAMEMESH_SLICE,itype,CNint,MaterialType1D,DATAIN,SLICEglo,rotMATglo,coorINT) ; 
             SLICE =  GeometryJoint(NAMEMESH_SLICE{itype},DATAIN) ;            
    end  
    SLICEglo{itype} = SLICE ;
    %
    IND_ELEM_TYPE = find(MaterialType1D == itype) ;  % Domains pertaining to type "itypes"
    CNintLOC = CNint(IND_ELEM_TYPE,:) ;
    
    nDOM = length(IND_ELEM_TYPE) ;   % Number of domains
    ndim = size(SLICE.COOR,2) ; % Spatial dimensions, 2 or 3
    if size(coorINT,2) < ndim & itype == 1
        coorINT = [coorINT, zeros(size(coorINT,1),1)] ;
    end
    % Distance between centroids
    DIST_CENTR = norm(SLICE.CENTRf2-SLICE.CENTRf1);
    nelem = size(SLICE.CN,1) ; % Number of elements per slice
    nnodeE = size(SLICE.CN,2) ; % Number of elements per slice
    nnode = size(SLICE.COOR,1) ; % Number of nodes per slice
    nnodeGLO = nDOM*nnode ;   % Total number of nodes
    nelemGLO = nDOM*nelem ;  % Total number of elements
    TypeElem = SLICE.TypeElement ; % Element type
    COORglo = zeros(ndim,nnodeGLO ) ; % Global coordinates
    CNglo = zeros(nelemGLO,nnodeE ) ; % Global coordinates
    MaterialTypeglo = zeros(nelemGLO,1) ; % Global coordinates
    
    rotMAT = zeros(ndim,ndim*size(CNintLOC,1)) ; % Global rotation matrix
    
    
    ifinCN = 0 ;ifinCOOR  = 0 ; ifin =0 ;
    if itype == 1
        iacumCN =   size(coorINT,1) ;
        
    end
    for e = 1:nDOM
        nodoINI = CNintLOC(e,1) ;
        nodoFIN = CNintLOC(e,2) ;
        xINI = coorINT(nodoINI,:) ;
        xFIN = coorINT(nodoFIN,:) ;
        DIST = norm(xFIN-xINI) ;
        if abs(DIST- DIST_CENTR)>1e-2*DIST_CENTR
            error(['The width of this elements is ',num2str(DIST_CENTR),' m' ])
        end
        % Rotation matrix
        r1 = (xFIN-xINI)'/DIST ;  % First vector
        r3 = [0,0,1]' ; % If not specified, the rotation is assumed to be around the x-axis
        r2 = cross(r3,r1) ;
        rotMATloc = [r1,r2,r3] ;
        iini =   ifin +1 ;
        ifin = iini + ndim -1;
        rotMAT(:,iini:ifin) = rotMATloc;
        %%%
        % Global coordinate matrix
        % -------------------------
        % Coordinates relative to centroid face 1
        COORrel  = zeros(size(SLICE.COOR')) ;
        for idim=1:ndim
            COORrel(idim,:) = SLICE.COOR(:,idim) - SLICE.CENTRf1(idim) ;
        end
        % Rotated coordinates
        COORrel = rotMATloc*COORrel ;
        % Translation
        for idim = 1:ndim
            COORrel(idim,:) = COORrel(idim,:) + xINI(idim) ;
        end
        % ------------
        iiniCO0R = ifinCOOR + 1;
        ifinCOOR = iiniCO0R + nnode -1;
        COORglo(:,iiniCO0R:ifinCOOR) = COORrel ;
        
        % Global connectivities
        iniCN = ifinCN + 1;
        ifinCN = iniCN + nelem -1 ;
        CNglo(iniCN:ifinCN,:) = SLICE.CN + iacumCN  ;
        
        MaterialTypeglo(iniCN:ifinCN) = SLICE.MaterialType+1;
        iacumCN  = iacumCN + nnode ;
    end
    
    rotMATglo{itype} = rotMAT ;
    
    
    if itype == 1
        COORprint = [coorINT; COORglo'] ;
        CNprint = {CNint,CNglo} ;
        TypeElementPRINT = {'Linear',SLICE.TypeElement} ;
        MaterialTypegloPRINT = {ones(size(CNint,1),1),MaterialTypeglo} ;
        NAMEMESH ={'1D',['3DTYPE ',num2str(itype)]} ;
        
    else
        COORprint = [COORprint; COORglo'] ;
        CNprint{itype+1} = CNglo ;
        TypeElementPRINT{itype+1} =  SLICE.TypeElement ;
        MaterialTypegloPRINT{itype+1} = MaterialTypeglo ;
        NAMEMESH{itype+1} = ['3DTYPE ',num2str(itype)] ;
    end
    
end

NAME_INPUT_DATA = [] ;
%GID PRINTING
[dummy1 NAMEMESH_STRUCTURE    dummy2]= fileparts(NAMEMESH_STRUCTURE) ;
NameFile_msh = ['GIDPOST/',NAMEMESH_STRUCTURE,'_',NAME_project,'.msh'] ;
IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,TypeElementPRINT,NAMEMESH);


% It proves necessary also to plot the direction of the local x-axis of each
% domain (1D plot). Remember that, by definition, this vector should point from face f1 to face f2
% -------------------------------------------------------------------------------------------------
NORMALS  =zeros(size(CNint,1),ndim) ;
%for e = 1:size(CNint,1)
NodeIni = CNint(:,1) ;
NodeFin = CNint(:,2) ;
COORini = coorINT(NodeIni,:) ;
COORfin = coorINT(NodeFin,:) ;
normCOOR = sqrt(sum((COORfin-COORini).^2,2)) ;
for idime = 1:ndim
    NORMALS(:,idime) = (COORfin(:,idime)-COORini(:,idime))./normCOOR ;
end
ElementsPrint = [] ;
for itype = 1:length(DATAIN.TYPE_SUBDOMAIN)
    switch DATAIN.TYPE_SUBDOMAIN{itype}
        case 'slice'
            ElementsPrint = [ElementsPrint; find(MaterialType1D==itype) ];
    end
end

%end
NameFile_res = ['GIDPOST/',NAMEMESH_STRUCTURE,'_',NAME_project,'.res'] ;
nelemE = [] ;
posgp = cell(2,1) ;
NODAL_VECTOR =[] ;
NODAL_SCALAR = [] ;
imesh = 1;
MESH(imesh).GAUSS_VECTOR(1).NAME = 'xlocal' ;
MESH(imesh).GAUSS_VECTOR(1).COMP = {'x1','x2','x3'}  ;
% All elements but joints

MESH(imesh).GAUSS_VECTOR(1).ELEMENTS =  ElementsPrint';
MESH(imesh).GAUSS_VECTOR(1).VAR = NORMALS(ElementsPrint,:)' ;
MESH(imesh).NAMEMESH = NAMEMESH{imesh};

GidResults2DFE_multi(NameFile_res,ndim,TypeElementPRINT,NODAL_VECTOR,NODAL_SCALAR,MESH,posgp);


%GidResults2DFE_multi(NameFile_res,ndim,nelemE,TypeElement,NODAL_VECTOR,NODAL_SCALAR,MESH,posgp);



