function [Vrb,Mintf,MESH,DATA,FluctuationFacesGLO] = FictInterfaceDISP_1D(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)

if nargin == 0
    load('tmp2.mat')
end
 % 1D interfaces (1D element)
 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/01_meta_1D.mlx
FluctuationFacesGLO = [] ; 
% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;

if isempty(DATAcommon.MESH_PARENT_DOMAIN_COARSE)
    error('yOU MUST DEFINE   DATAcommon.MESHcoarse')
else
    if  ischar(DATAcommon.MESH_PARENT_DOMAIN_COARSE)
        
        
         [MESHparent]= ReadMeshFileStr(DATAcommon.MESH_PARENT_DOMAIN_COARSE ,'READ_MATERIAL_COLUMN',0)  ;  
 % We only consider x-axis coordinates
 MESHparent.COOR = MESHparent.COOR(:,1) ; 
 
    else
        
        MESHparent.COOR = DATAcommon.MESH_PARENT_DOMAIN_COARSE.COOR ; 
        MESHparent.CN = DATAcommon.MESH_PARENT_DOMAIN_COARSE.CN ; 
    end
    
end
  
 
COORlin = MESHparent.COOR(MESHparent.CN,:) ;
CENTROID = MESH.GEOproperties.CENTROID ;
CENTROID = CENTROID(:,1) ; 
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORlin(:,idim) = COORlin(:,idim) - CENTROID(idim) ;
end

     

% 
% nnodeBND = 2;
% 
% nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
% ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
% DATAshape = ShapeFunCoefficients(COORlin,ORDER_POLYNOMIALS) ;
% DATAlocSHAPE.DATAshape  = DATAshape;
% xLIM = [] ;
% DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS; 
%[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% I think all these variables are dummy here... 
FESHAPE_coarse_elem_transf_coord.COOR =COORlin ; %
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORlin,1) ;
%FESHAPE_coarse_elem_transf_coord.DATA_ShapeFunctionFE =DATAlocSHAPE ;
%FESHAPE_coarse_elem_transf_coord.ShapeFunction = 'ShapeFunctionFE' ;


DATA.FESHAPE_coarse_elem_transf_coord = FESHAPE_coarse_elem_transf_coord;

% Vrb = [Vrb_face1,Vrb_face2]
ndim  = size(MESH.COOR,2) ; 
Vrb = zeros(ndim*length(BoundaryNodes),2) ; 
% FACE 
for iface = 1:2
NOD_LOC  =  MESH.BNDinterface(iface).NODES';
[~,IND_LOC_FACE,~] = intersect(BoundaryNodes,NOD_LOC) ;
DOFS_IND_LOC = small2large(IND_LOC_FACE,ndim) ; 
ModeLoc = zeros(length(DOFS_IND_LOC),1) ; 
ModeLoc(1:ndim:end) = 1; 
Vrb(DOFS_IND_LOC,iface) = ModeLoc ; 
end
 


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

 
% 
% PhiDEF_f = PhiDEF(MESH.faceDOFSall,:) ;
% PhiFLUC =   PprojDEF_operator(Vrb,Mintf,PhiDEF_f);

% ---------------------------------------------------------------
% REPLACE THESE MODES BY MODES WITH ZERO VALUE AT THE CORNERS
% ---------------------------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;

 


 
