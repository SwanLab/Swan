function [Vrb,DATA] = FictInterfaceModes_USER(MESH,DATAcommon,DATA,DATAoffline)
%--------------------------------------------------------------------------
% FUNCTION: FictInterfaceModes_USER
%
% PURPOSE:
%   Constructs interface displacement modes (`Vrb`) defined manually by the user
%   through a custom function specified in `DATAcommon.INFO_INTERFACE_MESH.FUN_TO_CONSTRUCT_MODES`.
%   These modes can be polynomial or user-defined shape functions and are used
%   to approximate the interfacial behavior of reduced-order or multiscale models.
%
% DESCRIPTION:
%   This function is invoked when the interface type is set to `'GIVEN_BY_USER'`
%   in the multiscale formulation (e.g., EIFEM). It evaluates the user-provided
%   function for generating shape functions (`Nshape`) on the boundary nodes of
%   the subdomain, referenced with respect to the centroid of the domain.
%   The resulting matrix of shape functions is used to construct the mode matrix `Vrb`.
%
%   The function also stores geometric and transformation information in the `DATA`
%   structure for later use in coarse-to-fine mapping, and exports the modes to
%   `.msh` and `.res` files for visualization in GiD.
%
% INPUTS:
%   - MESH        : Structure containing mesh data of the subdomain, including
%                   coordinates and boundary node lists (`faceNODESall`).
%   - DATAcommon  : High-level configuration structure. Must contain:
%                    - INFO_INTERFACE_MESH: structure with field
%                        .FUN_TO_CONSTRUCT_MODES: function handle for shape construction
%   - DATA        : Problem-specific structure. Enriched with fields for visualization
%                   and geometric transformation storage.
%   - DATAoffline : Structure with offline parameters (not used explicitly here).
%
% OUTPUTS:
%   - Vrb         : Matrix of size (ndof × nmodes), containing the user-defined
%                   displacement modes at the interface.
%   - DATA        : Updated structure with:
%                    • .FESHAPE_coarse_elem_transf_coord
%                    • .INFO_EDGES
%
% SIDE EFFECTS:
%   - Saves GiD mesh (`.msh`) and result (`.res`) files with the name:
%       'MODES/<DATA.NAME_BASE>DispIntfQuad.{msh,res}'
%
% DEPENDENCIES:
%   - RenumberConnectivities
%   - GidPostProcessModesDOML
%   - User-defined shape construction function: must have signature like
%       [Nshape,COORlin,INFO_EDGES] = FUN_TO_CONSTRUCT_MODES(INFO_INTERFACE_MESH, CENTROID, COORbnd)
%
% NOTES:
%   - Shape functions (`Nshape`) are assumed to be defined in local coordinates
%     with respect to the centroid of the domain.
%   - This function plays a key role in defining custom kinematic subspaces
%     for domain interfaces in ROM and multiscale frameworks.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Created: 2-April-2025, Barcelona
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------



% Modes fict. interfaces introduced by the user via a given function
% 2-APril-2025, UPC, CampusNord, BArcelona 

if nargin == 0
    load('tmp1.mat')
     
end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;

%% 
CENTROID = MESH.GEOproperties.CENTROID ; % ORIGIN Local reference system 
 
if ~isfield(DATAcommon,'INFO_INTERFACE_MESH')
    error('You must define the structure variable INFO_INTERFACE_MESH')
end
if ~isfield(DATAcommon.INFO_INTERFACE_MESH,'FUN_TO_CONSTRUCT_MODES')
    error('You must define the function for constructing the interface modes')
end 

disp('ADAPT THE CORRESPONDING N-SHAPE FUNCTION TO INCLUDE A GENERIC OUTPUT VARIABLE "INFO_EDGES", SEE  InputDataFunctions/LinearCubicPieceWise2Dmodes_EIFEM.m')
  

[Nshape,COORlin,INFO_EDGES] = feval(DATAcommon.INFO_INTERFACE_MESH.FUN_TO_CONSTRUCT_MODES,DATAcommon.INFO_INTERFACE_MESH,CENTROID,COORbnd) ; 

 FESHAPE_coarse_elem_transf_coord.COOR =COORlin ; %
 FESHAPE_coarse_elem_transf_coord.CN =1:size(COORlin,1) ;
DATA.FESHAPE_coarse_elem_transf_coord =   FESHAPE_coarse_elem_transf_coord;
DATA.INFO_EDGES = INFO_EDGES ;


%  
% COORlin = MESHparent.COOR(MESHparent.CN,:) ;
% for idim = 1:ndim
%     COORlin(:,idim) = COORlin(:,idim) - CENTROID(idim) ;
% end

     

% 
% nnodeBND = 4;
% 
% nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
% ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
% DATAshape = ShapeFunCoefficients(COORlin,ORDER_POLYNOMIALS) ;
% DATAlocSHAPE.DATAshape  = DATAshape;
% xLIM = [] ;
% DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
% [Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% Recall that this EIF element will play the role of "parent element" in
% our coares-scale formulation. Therefore, it is necessary to store the
% shape functions and the derivatives of the shape functions for purposes
% of transforming the coordinates when the physical domain does not
% coincide with the parent domain.
% 
% FESHAPE_coarse_elem_transf_coord.COOR =COORlin ; %
% FESHAPE_coarse_elem_transf_coord.CN =1:size(COORlin,1) ;
% FESHAPE_coarse_elem_transf_coord.DATA_ShapeFunctionFE =DATAlocSHAPE ;
% FESHAPE_coarse_elem_transf_coord.ShapeFunction = 'ShapeFunctionFE' ;






 ndim = size(CENTROID,2);

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
NameLoc =     'DispIntfQuad' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];
% ------------------------------------------
 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb ],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;

