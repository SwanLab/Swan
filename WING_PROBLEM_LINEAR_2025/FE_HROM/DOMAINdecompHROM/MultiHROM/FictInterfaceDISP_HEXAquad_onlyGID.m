function [Vrb,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_HEXAquad_onlyGID(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)
% Fictitious interface modes for hexahedra elements (27 nodes)
%  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/04_HIGHOR.mlx
% FictInterfaceDISP_HEXAquad turns out to be problematic,
% JAHO, 18-APRIL-2024
% ----------------------------------------
if nargin == 0
    load('tmp1.mat')
    
end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(DATAcommon.MESH_PARENT_DOMAIN_COARSE)
    error('yOU MUST DEFINE   DATAcommon.MESHcoarse')
end

 [MESHparent]= ReadMeshFileStr(DATAcommon.MESH_PARENT_DOMAIN_COARSE ,'READ_MATERIAL_COLUMN',0)  ; 


%MESHparent = DATAcommon.MESHcoarse ;
COORhexa = MESHparent.COOR(MESHparent.CN,:) ;
CENTROID = MESH.GEOproperties.CENTROID ;
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORhexa(:,idim) = COORhexa(:,idim) - CENTROID(idim) ;
end


%
DATAinp = [] ;
TypeElement = 'Hexahedra' ;

DATAoffline.InterpolationMethod = 'FEnodes'  ;

DATAcommon = DefaultField(DATAcommon,'ORDER_INVERSE_ELEMENT_inverse_mapping',4) ;

DATAoffline = DefaultField(DATAoffline,'ORDER_INVERSE_ELEMENT_inverse_mapping',DATAcommon.ORDER_INVERSE_ELEMENT_inverse_mapping) ;
DATAoffline = DefaultField(DATAoffline,'NPOINTS_ONE_DIRECTION_inverse_mapping',10) ;

%     DATAinp.ORDER_INVERSE_ELEMENT_inverse_mapping = DATAoffline.ORDER_INVERSE_ELEMENT_inverse_mapping;
%     DATAinp.NPOINTS_ONE_DIRECTION_inverse_mapping = DATAoffline.NPOINTS_ONE_DIRECTION_inverse_mapping;  ;

[Nshape,~,~] = ShapeFunDer_inversemapping(COORhexa',COORbnd,TypeElement,DATAoffline) ;
DATAlocSHAPE = []  ;
DATAlocSHAPE.DATAshape  = [];
%end




% These are the 9 shape functions corresponding to the 9 nodes of a
% quadratic element. However, the ninth shape function is the one
% corresponding to the "centroid" node, and therefore, its value at the
% boundary is identically zero ---thus, we have to get rid of the last
% column of Nshape
Nshape = Nshape(:,1:end-1) ;


%
% MESHLOC.COOR = COORhexa  ;
% MESHLOC.CN = [1:8];
% DATAloc.MESH = MESHLOC ;
% DATAloc.NumberOfGaussPointsPerElement = [2,2,2] ;
% DATAloc.MESH.TypeElement = 'Hexahedra' ;
% [~,DATA.weigthsGAUSSIAN,DATA.xGAUSSIAN] = GetMeshVariablesGAUSS(DATAloc) ;
%
%
%
%
% nnodeBND = length(CORNERS) ;
%
% nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
% ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
% DATAshape = ShapeFunCoefficients(COORhexa,ORDER_POLYNOMIALS) ;
% DATAlocSHAPE.DATAshape  = DATAshape;
% xLIM = [] ;
% DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
% [Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% Recall that this EIF element will play the role of "parent element" in
% our coares-scale formulation. Therefore, it is necessary to store the
% shape functions and the derivatives of the shape functions for purposes
% of transforming the coordinates when the physical domain does not
% coincide with the parent domain.

FESHAPE_coarse_elem_transf_coord.COOR =COORhexa ; %
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORhexa,1) ;
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
NameLoc =     'DispIntfhexa' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];
% ------------------------------------------

% -----------------
% FLUCTUATIONS
% ------------------
% ORTHOGONAL TO Vrb
% 
% if ~isempty(PhiDEF)
%     PhiDEF_f = PhiDEF(MESH.faceDOFSall,:) ;
%     PhiFLUC =   PprojDEF_operator(Vrb,Mintf,PhiDEF_f);
% else
%     PhiFLUC = [] ;
% end


% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb],[],NameFileMesh,NameFile_res,[],DATALOC) ;
INFOPLOTMESHBND.COOR = COORbnd;
INFOPLOTMESHBND.CN = CNbREN;
INFOPLOTMESHBND.TypeElement = MESH.TypeElementB;
INFOPLOTMESHBND.posgp = [];
INFOPLOTMESHBND.NameFileMesh = NameFileMesh;

 
