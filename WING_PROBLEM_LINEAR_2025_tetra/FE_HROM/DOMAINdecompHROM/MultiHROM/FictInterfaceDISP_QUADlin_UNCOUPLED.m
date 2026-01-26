function [Vrb,MESH,DATA] = FictInterfaceDISP_QUADlin_UNCOUPLED(MESH,DATAcommon,DATA,Mintf,DATAoffline)
% See
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.,
% Quadrilateral linear element, with uncoupled displacements in direction 1
% and 2 
% JAHO, 24-April-2024, Balmes 185, Barcelona
% -----------------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;

COORbnd = MESH.COOR(BoundaryNodes,:) ;

if isempty(DATAcommon.MESH_PARENT_DOMAIN_COARSE)
    error('yOU MUST DEFINE   DATAcommon.MESHcoarse')
end

DATA_LC.NameFileMeshDATA = DATAcommon.MESH_PARENT_DOMAIN_COARSE ; 
disp(['Reading parent domain mesh (recall that lines 1,2,3,4 must be identified, as well as materials (set to 1))'])
MESHparent = GeometryMesh(DATA_LC) ;

  
COORlin = MESHparent.COOR(MESHparent.CN,:) ;
CENTROID = MESH.GEOproperties.CENTROID ;
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORlin(:,idim) = COORlin(:,idim) - CENTROID(idim) ;
end

     


nnodeBND = 4;

nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
DATAshape = ShapeFunCoefficients(COORlin,ORDER_POLYNOMIALS) ;
DATAlocSHAPE.DATAshape  = DATAshape;
xLIM = [] ;
DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% Recall that this EIF element will play the role of "parent element" in
% our coares-scale formulation. Therefore, it is necessary to store the
% shape functions and the derivatives of the shape functions for purposes
% of transforming the coordinates when the physical domain does not
% coincide with the parent domain.

FESHAPE_coarse_elem_transf_coord.COOR =COORlin ; %
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORlin,1) ;
FESHAPE_coarse_elem_transf_coord.DATA_ShapeFunctionFE =DATAlocSHAPE ;
FESHAPE_coarse_elem_transf_coord.ShapeFunction = 'ShapeFunctionFE' ;

DATA.FESHAPE_coarse_elem_transf_coord = FESHAPE_coarse_elem_transf_coord;

nmodes = ndim*size(Nshape,2) ;
ndofsLOC = ndim*size(Nshape,1) ;
Vrb =zeros(ndofsLOC,nmodes) ;

for innode = 1:size(Nshape,2)
    for idim = 1:ndim
        imode = ndim*(innode-1)+idim ;
        Vrb(idim:ndim:end,imode) = Nshape(:,innode) ;
    end
end

 % The above are the mdoes of a standard quadrilateral element 
 % Now we have duplicate these modes, and for the first set of 8 modes,
 % set to zero the displacements of the faces in the direction 2, and for
 % the second set, the displacements in the direction 1 
 
 
 % Firstly, we have to determine the DOFs of each direction 

DATAcommon = DefaultField(DATAcommon,'Faces_Direction_Uncoupling',{[1,2],[3,4]});
ndim  = size(MESH.COOR,2) ; 
IND_loc_DOFs_faces_UNCOUP = cell(size(DATAcommon.Faces_Direction_Uncoupling)) ; 
IND_LOC_FACE = IND_loc_DOFs_faces_UNCOUP ; 
for IDIR = 1:length(DATAcommon.Faces_Direction_Uncoupling)     
    FACES_DIR = DATAcommon.Faces_Direction_Uncoupling{IDIR} ;
    NOD_LOC = cell(size(FACES_DIR))  ;
    for  iface = 1:length(FACES_DIR)
        NOD_LOC{iface} =  MESH.BNDinterface(FACES_DIR(iface)).NODES';
    end
    NOD_LOC = cell2mat(NOD_LOC') ;
    NOD_LOC = NOD_LOC(:);
    [~,IND_LOC_FACE{IDIR},~] = intersect(BoundaryNodes,NOD_LOC) ;
    IND_loc_DOFs_faces_UNCOUP{IDIR} = small2large(IND_LOC_FACE{IDIR},ndim);
end


%% FIRST SET OF MODES 
V_direction_1 = Vrb ; 
IPPOSITE = 2; 
V_direction_1(IND_loc_DOFs_faces_UNCOUP{IPPOSITE},:) = 0 ; 
V_direction_2 = Vrb ; 
IPPOSITE = 1; 
V_direction_2(IND_loc_DOFs_faces_UNCOUP{IPPOSITE},:) = 0 ; 
 

Vrb = [V_direction_1,V_direction_2] ; 




% PLOT  MODES
% -----------------
% Coordinates
% COORbnd = MESH.COOR(BoundaryNodes,:) ;
% Connectivities
CNb = cell(length(MESH.BNDinterface),1) ;
for iface  = 1:length(MESH.BNDinterface)
    CNb{iface} = MESH.BNDinterface(iface).CNb ;
end
CNb = cell2mat(CNb(:)) ;
CNbREN  =  RenumberConnectivities( CNb,1:length(BoundaryNodes) );


% ------------------------------------------------------------------------------------------------------------------------
NameLoc =     'DispIntfQuad' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];
% ------------------------------------------

% -----------------
% FLUCTU 

% ---------------------------------------------------------------
% REPLACE THESE MODES BY MODES WITH ZERO VALUE AT THE CORNERS
% ---------------------------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;
 