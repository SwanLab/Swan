function DATA_REFMESH =  PropertiesReferenceRVE(DATA_REFMESH,NameFileMesh,DOMAINVAR,COOR,...
    CN,CONNECTb,TypeElement,TypeElementB,MaterialType,posgp,density,DATAIN,NameFileMeshLOC_coarse)

if nargin == 0
    load('tmp1.mat')
end

idomREF = 1;
idomX = 1 ;
idomY = 1 ;
DATA_REFMESH.NameFileMesh = NameFileMesh ; % Name of the original mesh

% IDENTIFICATION OF FACE NODES
% -----------------------------
% TWO CASES --- CONTINUOUS STRUCTURES (FACES SHARING NODES)/ NON-CONTINUOUS
DATAIN  =DefaultField(DATAIN,'ContinuumStructuresWithoutCorners',0) ;
% Treat continuum structures without considering corners (23-May-2019)
% IMP-MAY19-A



if ~isempty(DOMAINVAR.NODES_CORNERS)
    %  CONTINUOUS STRUCTURES (FACES SHARING NODES)
    NODES_faces = DOMAINVAR.NODES_faces(idomX,idomY,:) ;
    NODES_faces = NODES_faces(:);
    DATA_REFMESH.NODES_faces = NODES_faces ; %
    
    % Continuous structures
    NODES_CORNERS = DOMAINVAR.NODES_CORNERS(idomX,idomY,:) ;
    NODES_CORNERS = NODES_CORNERS(:);
    DATA_REFMESH.NODES_CORNERS = NODES_CORNERS ; %
    %
    NODES_SIDES = DOMAINVAR.NODES_SIDES(idomX,idomY,:) ;
    NODES_SIDES = NODES_SIDES(:);
    DATA_REFMESH.NODES_SIDES = NODES_SIDES ; %
    
    % Boundary elements
    % ------------------
    CONNECTb_local = CONNECTb(idomX,idomY,:) ;
    CONNECTb_local = CONNECTb_local(:) ;
    DATA_REFMESH.CONNECTb = CONNECTb_local' ;
    
    if DATAIN.ContinuumStructuresWithoutCorners == 1
        
        DATA_REFMESH.NODES_CORNERS = [] ; 
        DATA_REFMESH.NODES_SIDES = [] ; 
         % CORNER 1 --> INT(4,1)
         % CORNER 2 --> INT(1,2)
         % CORNER 3 --> INT(2,3)
         % CORNER 4 --> INT(3,4)
         
         % CORNER 1 AND 2 --> FACE 1 
         % CORNER 3 AND 4 --> FACE 3 
         NODES_faces = NODES_SIDES ; 
         NODES_faces{1} = [ NODES_SIDES{1} ;NODES_CORNERS{1} ; NODES_CORNERS{2} ] ; 
         NODES_faces{3} = [ NODES_SIDES{3} ;NODES_CORNERS{4}  ; NODES_CORNERS{3}  ] ; 
         DATA_REFMESH.NODES_faces  = NODES_faces ; 
         
         %%% BOUNDARY CONNECTIVITIES HAVE TO BE ALSO MODIFIED 
         % -----------------------------------------------------
         
         FACESS = [2,4];
         IC={[2,3],[1,4]}  ;
         for ifaceLOC = 1:length(FACESS)
             iface = FACESS(ifaceLOC) ; 
             for iii  =1:length(IC{ifaceLOC})
                 ic = IC{ifaceLOC}(iii) ; 
                 [INDD,REMMM] = find(CONNECTb_local{iface} == NODES_CORNERS{ic }) ;
                 CONNECTb_local{iface}(INDD,:) = [] ;
             end
         end
          
          DATA_REFMESH.CONNECTb  = CONNECTb_local ; 
         
         
    end
    
    
else
    %  NON-CONTINUOUS
    NODES_faces = DOMAINVAR.NODES_faces(idomX,idomY,:) ;
    NODES_faces = NODES_faces(:);
    DATA_REFMESH.NODES_faces = NODES_faces ; %
    
    % Boundary elements
% ------------------
CONNECTb_local = CONNECTb(idomX,idomY,:) ;
CONNECTb_local = CONNECTb_local(:) ;
DATA_REFMESH.CONNECTb = CONNECTb_local' ;
    
end


% Coordinates
ListNodes = DOMAINVAR.ListNodesDom{idomX,idomY} ;
DATA_REFMESH.COOR = COOR(ListNodes,:) ;

DATA_REFMESH.NameFileMeshLOC_coarse = NameFileMeshLOC_coarse ; % Name COARSE MESH


%%%%%%%%%%%%%%%%%%%%
if strcmp(NameFileMesh,NameFileMeshLOC_coarse) ==0  && ~isempty(NameFileMeshLOC_coarse)
    % Determining new connectivities
    [DATA_REFMESH.CONNECTb_coarse DATA_REFMESH.TypeElementB_coarse]= ...
        MatchingFineCoarseMeshRVE(NameFileMeshLOC_coarse,DATA_REFMESH.COOR ) ;
else
    DATA_REFMESH.CONNECTb_coarse = [] ;
    DATA_REFMESH.TypeElementB_coarse = [] ;
end



% Connectivities
ELEMS =  DOMAINVAR.ListElements{idomREF} ;
DATA_REFMESH.CN = CN(ELEMS,:) ;
% How to be sure that all entries of CN are consecutive  ?
ALLNODES = unique(DATA_REFMESH.CN) ;
nnodeREF = size(DATA_REFMESH.COOR,1 ) ;
ALLNODES_REF = 1:nnodeREF;
if norm(ALLNODES-ALLNODES_REF') ~= 0
    % Renumbering is required
    DATA_REFMESH.CN = RenumberConnectivities(DATA_REFMESH.CN,ALLNODES_REF) ;
end

% Other mesh information
DATA_REFMESH.TypeElement = TypeElement ;
DATA_REFMESH.TypeElementB = TypeElementB ;
DATA_REFMESH.posgp = posgp ;
% Material type
DATA_REFMESH.MaterialType = MaterialType(ELEMS) ;
DATA_REFMESH.density = density(ELEMS) ;

% It only remains to compute the global matrix of boundary shape
% functions. This global matrices will allow one to fastly compute the
% nodal equivalent forces caused by distributed loads over each of the
% surfaces of the slice. The idea is that the nodal forces in the direction i,
% denoted by F_idim (as many entries as nodes), may be computed as

%  F_idim = N_left{1}'*W{1}*Tgauss_idim{1}  +  N_left{2}'*W{2}*Tgauss_idim{2} + ...N_left{isurf}'*W{isurf}*Tgauss_idim{isurf}
%  Here, Tgauss_idim{isurf} is the vector of traction forces at the
%  Gauss points of boundary elements of surface isurf
% In turn, if we only know the value of the traction forces at the
% corresponding nodes of the surface, then we can make
% Tgauss_idim{isurf} = N_right{isurf}*Tnodes_idim{isurf}

% In fact, if we can make N_left{isurf} = N_right{isurf} = Nbnd{isurf}. To this end,
% Tnodes_idim{isurf} should be defined for all nodes of the domain. We
% shall use this latter option
nfaces = length(DATA_REFMESH.CONNECTb) ;
Nbnd = cell(1,nfaces) ;
Wbnd = cell(1,nfaces) ;

for iface = 1:nfaces
    
    % Compute boundary shape functions of surface "iface"
    CNb = DATA_REFMESH.CONNECTb{iface} ;
    [ NelemB,wSTb ] = ComputeNelemBoundALL(DATA_REFMESH.COOR ,CNb,TypeElementB) ;
    % NelemB = nelemB*ngaus x nnode  --> Shape functions
    % wSTb --> Gauss weights (including Jacobians )
    nnode = size(DATA_REFMESH.COOR,1);  % Number of total nodes
    nelemB = size(CNb,1);  % Number of boundary elements (total)
    nnodeEb = size(CNb,2) ; % Number of nodes per boundary element
    ndim = size(DATA_REFMESH.COOR,2); % Number of spatial dimensions
    ngaus = size(NelemB,1)/nelemB ;
    Nbnd{iface} = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngaus,CNb,nnode) ;
    % Nbnd{iface} is a nelemB*ngaus x nnode  matrix. It relates the values
    % of nodal variables defined at all nodes to values at the Gauss points
    % of the boundary nodes of CNb
    ndimLOC = 1;
    Wbnd{iface} = CompWeightDiag(wSTb,ndimLOC)  ; % Diagonal weight matrix
end

DATA_REFMESH.Nbnd = Nbnd ;
DATA_REFMESH.Wbnd = Wbnd ;


% ROTATIONS ANGLE INTERFACES ----RELATIVE TO FACE 1 
%%%% Lastly, we compute the  rigid body modes of the interfaces
%  These variables are given with respect to the local
% system attatched to each face. This reference system coincides with the
% reference system of the domain itself for face 1.  

% RECALL THAT: FACE1 --> XMIN,  FACE3 --> XMAX
%              FACE2 --> YMIN,  FACE4 ---> XMAX    

 nfaces=  length(DATA_REFMESH.NODES_faces) ; 

 DATA_REFMESH.RotationMatrixFace  = cell(1,nfaces) ; 
 
 
% FACE 3 IS ROTATED WITH RESPECT TO FACE 1 (AROUND THE Y AXIS)
RotF3 = [] ;
if ~isempty(DATAIN.angDOM)
    angDOM = DATAIN.angDOM ;
    RotF3 = [     
           cos(angDOM)  0  sin(angDOM)  
             0          1     0
        -sin(angDOM)    0   cos(angDOM) ] ;
    
end
DATA_REFMESH.RotationMatrixFace{3} =RotF3; 


for iface = 1:nfaces
    [R,GGG,CentroidFA,Mst] = GeometricDataFaceRVE(iface,DATA_REFMESH,TypeElementB) ;
    DATA_REFMESH.CENTRf{iface} = CentroidFA ;
    DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{iface} = GGG ;
    DATA_REFMESH.GeometricMassMatrixInterface{iface} = Mst ;
    DATA_REFMESH.RigidBodyModesInterface{iface} = R ;
end

%%% CONTINUOUS DOMAINS. MODES OF CORNERS, AND SIDES 1 AND 2 (= 3 AND 4)
% ----------------------------------------------------------------------
DATA_REFMESH = DefaultField(DATA_REFMESH,'NODES_CORNERS',[]) ;
if ~isempty(DATA_REFMESH.NODES_CORNERS)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DATA_REFMESH = CornersSidesProperties_PLATES(DATA_REFMESH) ;
    
    
    
else
    
    DATA_REFMESH.ModesSides = [] ;
    DATA_REFMESH.ModesCorners = [] ;
    
    
    
end





DATAIN = DefaultField(DATAIN,'MASS_MATRIX_FORMULATION_REACTIONS',1) ;

if  DATAIN.MASS_MATRIX_FORMULATION_REACTIONS == 1
    %load(DATAIN.NAME_WS_MODES,'Nst') ;
    DATA_REFMESH.M =  DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ; % speye(size(DATA_REFMESH.Nst,2));% DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ;
    
else
    
    DATA_REFMESH.M =   speye(ndim*nnode );% DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ;
    
end

end

