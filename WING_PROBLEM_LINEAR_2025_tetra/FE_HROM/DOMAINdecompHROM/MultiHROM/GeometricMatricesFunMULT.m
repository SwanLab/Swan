function [OPERFE,MESH] =  GeometricMatricesFunMULT(MESH,MATPRO)
%--------------------------------------------------------------------------
% FUNCTION: GeometricMatricesFunMULT
%
% PURPOSE:
%   Constructs and stores geometric and finite element operators for a
%   given subdomain mesh. These operators are precomputable, meaning they
%   do not depend on the specific displacement solution and are therefore
%   reusable across simulations or in model reduction contexts (e.g., EIFEM).
%
%   This function is adapted for multiscale (RVE-based) settings, including:
%   - Strain-displacement matrix `Bst`
%   - Shape function matrix `Nst`
%   - Mass matrix `M`
%   - Gauss weights and integration points
%   - Boundary shape functions and geometric mass matrices
%   - Domain centroid and rotational inertias
%   - Rigid body modes (domain and faces)
%
% USAGE:
%   [OPERFE, MESH] = GeometricMatricesFunMULT(MESH, MATPRO)
%
% INPUT:
%   - MESH: structure with domain geometry, connectivities, and element types
%   - MATPRO: material properties (only density is used here for mass matrix)
%
% OUTPUT:
%   - OPERFE: structure with assembled FE operators and matrices
%   - MESH  : updated mesh with geometric properties and local coordinates
%
% WHAT THIS FUNCTION DOES:
%   1. Computes B-matrix for small or large strains (depending on flags)
%   2. Computes shape function matrix `Nst` at Gauss points for RHS/mass
%   3. Assembles consistent mass matrix `M` using material density
%   4. Computes rigid body modes and centroid of the domain
%   5. Translates coordinates so the domain centroid is the origin
%   6. Assembles boundary shape function matrices:
%      - Right operator: maps nodal values to Gauss points
%      - Left operator: projects Gauss-point values back to nodal level
%   7. Computes weighted products for each boundary face defined in GID
%   8. For each face:
%      - Computes geometric mass matrix
%      - Computes local centroids and relative coordinates
%      - Stores normals and tangents at Gauss points
%
% REMARKS:
%   - The function supports multiscale models where geometric metrics of
%     both the domain and interface faces must be explicitly precomputed.
%   - Rigid body modes are computed via mass-weighted basis construction.
%   - Data in `OPERFE` will be used in stiffness assembly and projection.
%
% SEE ALSO:
%   - BstLargeStrains, BstSmallStrains
%   - ComputeNelemALL, AssemblyNGlobal
%   - CentroidGeometricMassMatrixNEW
%   - ConstructBasisRigidBody_MASSM
%   - NormalsAndTangentVectorsGaussPoints
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC - CIMNE, Barcelona
%   Date: 09-Feb-2023
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------




% Matrix depending on the finite element mesh (i.e., that can be pre-computed)
% Copy of GeometricMatricesFun.m. Adapted to   MULTISCALE problems (RVE)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx

% For instance, for a problem with linear quadrilateral elements (body), 100 body
% elements, 121 nodes, 40 boundary elements (2 Gauss points per element),
% and 4 faces: 

% MESH = 
% 
%   struct with fields:
% 
%                            COOR: [121×2 double]
%                              CN: [100×4 double]
%                    MaterialType: [100×1 double]
%                     TypeElement: 'Quadrilateral'
%                     NODES_FACES: {[11×1 double]  [11×1 double]  [11×1 double]  [11×1 double]}
%                             CNb: [40×2 double]
%                    TypeElementB: 'Linear'
%      NormalBoundaryElementsFace: {[2×10 double]  [2×10 double]  [2×10 double]  [2×10 double]}
%     TangentBoundaryElementsFace: {[2×10 double]  [2×10 double]  [2×10 double]  [2×10 double]}
%                            DATA: []
%                     posgp_given: []
%       Indexes_faces_bnd_element: {[1 2 3 4 5 6 7 8 9 10]  [11 12 13 14 15 16 17 18 19 20]  [21 22 23 24 25 26 27 28 29 30]  [31 32 33 34 35 36 37 38 39 40]}
%                   GEOproperties: [1×1 struct]
%                PROPERTIES_FACES: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}

% OPERFE = 
%
% 
%                    Bst: [1600×242 double]
%                   wSTs: [400×1 double]
%           ngaus_STRESS: 4
%                  posgp: [2×4 double]
%               wSTs_RHS: [400×1 double]
%              posgp_RHS: [2×4 double]
%              ngaus_RHS: 4
%                    Nst: [800×242 double]
%                   wSTb: [80×1 double]
%              NstB_left: [80×121 double]
%    NstT_W_N_boundaries: {[121×20 double]  [121×20 double]  [121×20 double]  [121×20 double]}
% See more info in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% --------------------------------------------------------
% JAHO - 9-Feb-2023  
if nargin == 0
    load('tmp.mat')
end
DATA = MESH.DATA;

nstrain = DATA.MESH.nstrain ; 

%1 )  Bmatrix. Matrix such that DeformationGradient = Identity + Bst_F*d
% wSTs --- Gauss weights x Jacobian

% This is a variable to disable the use of the deformation gradient in
% Small strain kinematics
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0); 
DATA = DefaultField(DATA,'NO_USE_Deformation_gradient_in_Small_Strains',0); 

if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1 && DATA.SMALL_STRAIN_KINEMATICS ==1
    [OPERFE.Bst,OPERFE.wSTs,OPERFE.ngaus_STRESS,IDENTITY_F,OPERFE.posgp] = BstSmallStrains(MESH,nstrain) ;
else
    [OPERFE.Bst,OPERFE.wSTs,OPERFE.ngaus_STRESS,IDENTITY_F,OPERFE.posgp] = BstLargeStrains(MESH,nstrain) ;
end
MESH.DATA  = [] ;
 

% Stacked N-matrix. A matrix such that  u(x)  = Nst*d, where u(x) is the displacement at a given Gauss Point
nnode = size(MESH.COOR,1); ndim = size(MESH.COOR,2); nelem = size(MESH.CN,1); nnodeE = size(MESH.CN,2) ;
% Number of Gauss points for integrating RHS and mass MATRIX may be different
% from number of Gauss points used in integrating stiffness matrix

MESH = DefaultField(MESH,'posgp_given',[]) ;
if isempty(MESH.posgp_given)
    TypeIntegrand = 'RHS';
else
    TypeIntegrand = {MESH.posgp_given,MESH.weights_given} ;
end

[~,~,shapef_RHS,~] = ComputeElementShapeFun(MESH.TypeElement,nnodeE,TypeIntegrand) ;

[  OPERFE.wSTs_RHS,~, OPERFE.posgp_RHS] = ComputeW_RHS(MESH.COOR,MESH.CN,MESH.TypeElement,ndim,TypeIntegrand)  ;

[ Nshapeelem,~  ] = ComputeNelemALL(MESH.TypeElement,nnodeE,ndim,nelem,TypeIntegrand) ;
OPERFE.ngaus_RHS = size(OPERFE.posgp_RHS,2) ;
OPERFE.Nst = AssemblyNGlobal(Nshapeelem,nelem,nnodeE,ndim,OPERFE.ngaus_RHS,MESH.CN,nnode) ;
%wDIAG_RHS = CompWeightDiag(wSTs_RHS,ndim)  ;

% 5.  MASS MATRIX
% -----------------------
densGLO = repmat(MATPRO.dens',DATA.MESH.ngaus_RHS,1) ;
densGLO = densGLO(:) ;
densGLO_W = CompWeightDiag(densGLO.*OPERFE.wSTs_RHS,ndim) ;
NstW = densGLO_W*OPERFE.Nst ;
OPERFE.M = NstW'*OPERFE.Nst ;


% COMPUTING THE CENTROID OF THE DOMAIN, AS WELL AS THE ROTATIONAL INERTIAS.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
DATALOC = [] ;
[BasisUrb,VOLUME,Rbar,CENTROID,INERTIA ]= ConstructBasisRigidBody_MASSM(MESH.COOR,OPERFE.Nst,OPERFE.wSTs_RHS,DATALOC) ; %

GEOproperties.RIGID_BODY_MODES = BasisUrb ;
GEOproperties.CENTROID = CENTROID ;
GEOproperties.INERTIA = INERTIA ;
GEOproperties.VOLUME = VOLUME ;

%   COORDINATES with respect to centroid
% ---------------------------
for idim = 1:size(MESH.COOR,2)
    MESH.COOR(:,idim) = MESH.COOR(:,idim) - CENTROID(idim) ; 
end

% BOUNDARY OPERATORS
CNb = MESH.CNb ; 
MESH.CNb = [] ; 
Indexes_faces_bnd_element = cell(size(CNb)) ; 
iacum = 1; 
for ifaces  = 1:length(CNb) ; 
    MESH.CNb = [MESH.CNb;  CNb{ifaces}] ; 
    Indexes_faces_bnd_element{ifaces} = iacum:(iacum+size(CNb{ifaces},1)-1) ; 
    iacum = iacum + size(CNb{ifaces},1) ; 
end
MESH.Indexes_faces_bnd_element = Indexes_faces_bnd_element; 


% --------------------BOUNDARY MATRICES

%disp('Computing shape function matrices for all boundary elements... (Nst)')
[ NelemB ,OPERFE.wSTb ] = ComputeNelemBoundALL(MESH.COOR,MESH.CNb,MESH.TypeElementB) ;
% -------------------------------------------------------------------------
% Assembly of matrix NstB. This is a matrix relating the nodal forces at
% the boundary elements with the nodal forces at the  Gauss points of the
% boundary elements
nelemB = size(MESH.CNb,1);  % Number of boundary elements (total)
nnodeEb = size(MESH.CNb,2) ; % Number of nodes per boundary element
ngausB = size(NelemB,1)/nelemB ;
%disp('Assembly of NstB...')
NstB_right = AssemblyNboundRIGHT(NelemB,nelemB,nnodeEb,ndim-1,ngausB,nnode) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Left-operator NstBw'*NstB*Tnod
%%%
% This is a matrix relating the forces/displacements at all nodes of the
% discretization, with the values of the forces at the Gauss points of the
% boundary elements
disp('Assembly of NstBw...')
OPERFE.NstB_left = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngausB,MESH.CNb,nnode) ;
% Diagonal matrix with weights
ndimLOC = 1;
wDIAGb = CompWeightDiag(OPERFE.wSTb,ndimLOC)  ;
NstB_leftW = wDIAGb*OPERFE.NstB_left ;
% f = unique(MESH.CNb(:)) ; 
% Mbound = NstB_left'*NstB_leftW ; 
% Mff = Mbound(f,f) ; 

% We are only interested in the contribution of such matrices at the faces
% defined in GID pre-process. Accordingly, we make
NstT_W_N_boundaries = cell(size(MESH.NODES_FACES)) ;

for iface = 1:length(NstT_W_N_boundaries)
    indlocal=  MESH.Indexes_faces_bnd_element{iface} ;
    indlocal_columns = small2large(indlocal,nnodeEb)  ;
    indlocal_rows =  small2large(indlocal,ngausB)  ;
    NstT_W_N_boundaries{iface} = NstB_leftW(indlocal_rows,:)'*NstB_right(indlocal_rows,indlocal_columns) ;
    
end
OPERFE.NstT_W_N_boundaries = NstT_W_N_boundaries; 


% Rigid body modes and area of each of the faces defined by the user in GID
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FACES = cell(size(MESH.Indexes_faces_bnd_element)) ;
for iface = 1:length(FACES)
    indELEMB = MESH.Indexes_faces_bnd_element{iface} ;
    CONNECTb_iface = MESH.CNb(indELEMB,:) ;
    [CentroidFA,AREA,Mst,Nst_face,wSTb_face] =CentroidGeometricMassMatrixNEW(MESH.COOR,[],CONNECTb_iface,MESH.TypeElementB) ;
    nodes  = unique(CONNECTb_iface(:)) ;
    COOR_FACE = MESH.COOR(nodes,:) ; % Coordinates of this face
    COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
    %  BasisRrb = ConstructBasisRigidBody(COORrelA) ;
    FACES{iface}.AREA = AREA ;
    %   FACES{iface}.RIGID_BODY_MODES = BasisRrb ;
    FACES{iface}.COORrelA_global = COORrelA ;
    FACES{iface}.CENTROID = CentroidFA ;
    FACES{iface}.GeometricMassMatrix = Mst ;
        FACES{iface}.Nst  = Nst_face ;
        FACES{iface}.wST  = wSTb_face ;

    
    
    [unitNORMALS,tTANGiniST] = NormalsAndTangentVectorsGaussPoints(MESH,CONNECTb_iface) ;
    
    FACES{iface}.UnitNormalAtGaussPoint = unitNORMALS ;
    FACES{iface}.UnitTangentAtGaussPoint = tTANGiniST ;
    
    
end
MESH.GEOproperties = GEOproperties; 
MESH.PROPERTIES_FACES = FACES; 
