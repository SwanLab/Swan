function DATA_REFMESH =  PropertiesReferenceDomain(DATA_REFMESH,NameFileMesh,DOMAINVAR,COOR,faceNODES_1,faceNODES_2,...
    CN,CONNECTb,TypeElement,TypeElementB,MaterialType,posgp,density,DATAIN,NameFileMeshLOC_coarse)

if nargin == 0
    load('tmp4.mat')
end




idomREF = 1;
DATA_REFMESH.NameFileMesh = NameFileMesh ; % Name of the original mesh


DATA_REFMESH.NODES_faces12 = DOMAINVAR.NODES_faces12(idomREF,:) ; % Nodes faces 1 and 2 (sorted)
% Length reference domain
%x_1 = COOR(faceNODES_1(1),1) ;
%x_2 =  COOR(faceNODES_2(1),1) ;
%DATA_REFMESH.LENGTH = x_2-x_1 ;
% Coordinates
DATA_REFMESH.COOR = COOR(DOMAINVAR.ListNodesDom{idomREF},:) ;

DATA_REFMESH.NameFileMeshLOC_coarse = NameFileMeshLOC_coarse ; % Name COARSE MESH


%%%%%%%%%%%%%%%%%%%%
if strcmp(NameFileMesh,NameFileMeshLOC_coarse) ==0  && ~isempty(NameFileMeshLOC_coarse)
    % Determining new connectivities
    [DATA_REFMESH.CONNECTb_coarse DATA_REFMESH.TypeElementB_coarse]= ...
        MatchingFineCoarseMesh(NameFileMeshLOC_coarse,DATA_REFMESH.COOR ) ;
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
if norm(ALLNODES(:)-ALLNODES_REF(:)) ~= 0
    % Renumbering is required
    DATA_REFMESH.CN = RenumberConnectivities(DATA_REFMESH.CN,ALLNODES_REF) ;
end
% Boundary elements
% ------------------
DATA_REFMESH.CONNECTb = CONNECTb(idomREF,:) ;
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


%%%% Lastly, we compute the  rigid body modes of the interfaces
% JAHO, 26-Nov-2018. These variables are given with respect to the local
% system attatched to each face. This reference system coincides with the
% reference system of the domain itself for face 1. For face 2, on curved
% domains, the relationship between the domain reference system and the
% face reference system is given by the matrix: RotF2 
RotF2 = [] ; 
DATAIN = DefaultField(DATAIN,'ISTWIST_ANGLE',0) ; 
if ~isempty(DATAIN.angDOM)
    angDOM = DATAIN.angDOM ; 
    if size(COOR,2) == 3
        if DATAIN.ISTWIST_ANGLE == 0
        RotF2 = [cos(angDOM)   sin(angDOM)  0
            -sin(angDOM) cos(angDOM)   0
            0            0        1 ] ;
        else
             RotF2 = [1 0 0
                 0 cos(angDOM)   -sin(angDOM)  
           0    sin(angDOM) cos(angDOM)   ] ;
        end
    else
        RotF2 = [cos(angDOM) sin(angDOM)
            -sin(angDOM)       cos(angDOM)] ;
    end
end

% ROTATION MATRIX of face 2 with respect to face 1 (the latter is assumed to be parallel to the global coordinates)
% ----------------------------------


%%%% FACE 1
iface=1 ;
DATA_REFMESH.RotationMatrixFace{1} = [] ; 
[R,GGG,CentroidFA,Mst,DOFA,Rwarping] = GeometricDataFace(iface,DATA_REFMESH,TypeElementB) ;
DATA_REFMESH.CENTRf1 = CentroidFA ;
DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{iface} = GGG ;
DATA_REFMESH.GeometricMassMatrixInterface{iface} = Mst ;
DATA_REFMESH.RigidBodyModesInterface{iface} = R ;
DATA_REFMESH.WarpingModesInterface{iface} = Rwarping ;

% FACE 2
% ---------
iface=2 ;
 DATA_REFMESH.RotationMatrixFace{2} =RotF2 ; 
     
[R,GGG,CentroidFB,Mst,DOFB] = GeometricDataFace(iface,DATA_REFMESH,TypeElementB) ;
DATA_REFMESH.CENTRf2 = CentroidFB ;
DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{iface} = GGG ;
DATA_REFMESH.GeometricMassMatrixInterface{iface} = Mst ;
DATA_REFMESH.RigidBodyModesInterface{iface} = R ;
DATA_REFMESH.WarpingModesInterface{iface} = Rwarping ;


DATA_REFMESH.LENGTH = norm(CentroidFB-CentroidFA) ;

% Mass matrix
% -----------
%DATA_REFMESH.M = DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ;

DATAIN = DefaultField(DATAIN,'MASS_MATRIX_FORMULATION_REACTIONS',1) ;

if  DATAIN.MASS_MATRIX_FORMULATION_REACTIONS == 1
  %  load(DATAIN.NAME_WS_MODES,'Nst') ;
    DATA_REFMESH.M =  DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ; % speye(size(DATA_REFMESH.Nst,2));% DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ;
  
    if norm(DATA_REFMESH.M(DOFA,DOFB),'fro')/norm(DATA_REFMESH.M(DOFA,DOFA),'fro') >1e-6
        warning('Option not valid for 1-element width slices')
        disp(['set DATAIN.MASS_MATRIX_FORMULATION_REACTIONS = 0'])
        error(' ')
    end
    
    % error('This option gives rise to problems in the reconstruction of rigid motions')
else
    
    DATA_REFMESH.M =   speye(ndim*nnode );% DATA_REFMESH.Nst'*(DATA_REFMESH.WdiagRHS*DATA_REFMESH.Nst) ;
    
end
