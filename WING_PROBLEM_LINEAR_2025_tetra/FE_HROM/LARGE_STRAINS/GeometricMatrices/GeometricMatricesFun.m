function [Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,...
    posgp,shapef_RHS] = ...
    GeometricMatricesFun(MESH,nstrain)
%%
% =========================================================================
% GeometricMatricesFun — Precompute geometric/shape operators & quadrature
% =========================================================================
% PURPOSE
%   Build all mesh-dependent operators that can be precomputed before the FE
%   solve: deformation-gradient operator Bst_F, global stacked shape matrix
%   Nst, quadrature weights/points for stiffness/RHS, boundary transfer
%   operators for tractions, and basic geometric properties (centroid,
%   volume, inertia, face normals/tangents, rigid-body modes).
%
% SIGNATURE
%   [Bst_F, wSTs, Nst, wSTs_RHS, NstT_W_N_boundaries, ngaus_RHS, ...
%    GEOproperties, ngaus_STRESS, IDENTITY_F, posgp, shapef_RHS] = ...
%       GeometricMatricesFun(MESH, nstrain)
%
% INPUTS
%   MESH        (struct)  Mesh data: COOR, CN, TypeElement, boundary meshes
%                         (CNb, TypeElementB), optional fields:
%                           • DATA.SMALL_STRAIN_KINEMATICS (logical)
%                           • DATA.NO_USE_Deformation_gradient_in_Small_Strains (logical)
%                           • posgp_given / weights_given (override RHS quadrature)
%   nstrain     (int)     Length of strain-like vector used for stresses.
%
% OUTPUTS
%   Bst_F           : Operator s.t. F ≈ I + Bst_F * d   (or small-strain B)
%   wSTs            : Gauss weights (× Jacobians) for stiffness/stress loops
%   Nst             : Global stacked shape matrix for nodal fields at RHS quad
%   wSTs_RHS        : Gauss weights for RHS/mass integration
%   NstT_W_N_boundaries : Cell array with boundary transfer operators
%   ngaus_RHS       : # of RHS Gauss points per element set
%   GEOproperties   : Struct with CENTROID, INERTIA, VOLUME, RIGID_BODY_MODES,
%                     boundary weights (wSTb), left operator NstB_left, and FACES info
%   ngaus_STRESS    : # of stress/stiffness Gauss points per element set
%   IDENTITY_F      : Identity part used with Bst_F in large/small strain forms
%   posgp           : Stiffness/stress quadrature points (global list)
%   shapef_RHS      : Element-level shape functions at RHS quadrature
%
% WHAT THE FUNCTION DOES
%   1) Chooses deformation-gradient operator:
%        • If SMALL_STRAIN_KINEMATICS==1 AND NO_USE_Deformation_gradient_in_Small_Strains==1
%          → BstSmallStrains(MESH,nstrain)
%        • Else → BstLargeStrains(MESH,nstrain)
%      Returns Bst_F, wSTs, ngaus_STRESS, IDENTITY_F, posgp.
%   2) Builds RHS/mass integration data:
%        • If MESH.posgp_given is empty → default RHS rule; otherwise use the provided
%          points/weights.
%        • Computes shapef_RHS, wSTs_RHS, posgp_RHS, and assembles global Nst.
%   3) Computes global geometric properties via ConstructBasisRigidBody_MASSM:
%        • Rigid-body modes (mass-orthonormal), total VOLUME, CENTROID, INERTIA tensor.
%   4) Assembles boundary operators:
%        • Compute element-level boundary shapes/weights (NelemB, wSTb).
%        • Build left/right boundary operators and weight them (NstB_leftW).
%        • For each user face, extract sub-operators → NstT_W_N_boundaries{iface}.
%   5) Per-face geometric data:
%        • CENTROID, AREA, geometric mass matrix, unit normals and tangents
%          at Gauss points (NormalsAndTangentVectorsGaussPoints).
%
% KEY DEPENDENCIES
%   BstSmallStrains / BstLargeStrains, ComputeElementShapeFun, ComputeW_RHS,
%   ComputeNelemALL, AssemblyNGlobal, ConstructBasisRigidBody_MASSM,
%   ComputeNelemBoundALL, AssemblyNboundRIGHT/LEFT, CompWeightDiag,
%   CentroidGeometricMassMatrixNEW, NormalsAndTangentVectorsGaussPoints,
%   small2large.
%
% PRACTICAL NOTES
%   • Stiffness/stress and RHS quadrature can differ (wSTs vs wSTs_RHS).
%   • Boundary transfer operators map Gauss-point tractions on selected faces
%     to the global nodal force vector (consistent loading).
%   • posgp_given / weights_given let you enforce custom RHS quadrature.
%   • GEOproperties.FACES{iface} carries AREA, CENTROID, normals/tangents, and
%     per-face geometric mass; useful for postprocessing and hydro loads.
%
% USAGE EXAMPLE
%   [Bst_F,wSTs,Nst,wSTs_RHS,Nb,ngRHS,GEO,ngST,I,posgp,shRHS] = GeometricMatricesFun(MESH, nstrain);
%   % Use Bst_F & wSTs for stresses/stiffness; Nst & wSTs_RHS for mass/RHS.
%
% .mlx references: (none)
%
% Dates/places:
%   Original note: JAHO — 22-Nov-2020 (Sunday)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

if nargin == 0
    load('tmp.mat')
end
DATA = MESH.DATA;


%1 )  Bmatrix. Matrix such that DeformationGradient = Identity + Bst_F*d
% wSTs --- Gauss weights x Jacobian

% This is a variable to disable the use of the deformation gradient in
% Small strain kinematics
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0); 
DATA = DefaultField(DATA,'NO_USE_Deformation_gradient_in_Small_Strains',0); 

if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1 && DATA.SMALL_STRAIN_KINEMATICS ==1
    [Bst_F,wSTs,ngaus_STRESS,IDENTITY_F,posgp] = BstSmallStrains(MESH,nstrain) ;
else
    [Bst_F,wSTs,ngaus_STRESS,IDENTITY_F,posgp] = BstLargeStrains(MESH,nstrain) ;
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

[  wSTs_RHS,~, posgp_RHS] = ComputeW_RHS(MESH.COOR,MESH.CN,MESH.TypeElement,ndim,TypeIntegrand)  ;

[ Nshapeelem,~  ] = ComputeNelemALL(MESH.TypeElement,nnodeE,ndim,nelem,TypeIntegrand) ;
ngaus_RHS = size(posgp_RHS,2) ;
Nst = AssemblyNGlobal(Nshapeelem,nelem,nnodeE,ndim,ngaus_RHS,MESH.CN,nnode) ;
%wDIAG_RHS = CompWeightDiag(wSTs_RHS,ndim)  ;

% COMPUTING THE CENTROID OF THE DOMAIN, AS WELL AS THE ROTATIONAL INERTIAS.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
DATALOC = [] ;
[BasisUrb,VOLUME,Rbar,CENTROID,INERTIA ]= ConstructBasisRigidBody_MASSM(MESH.COOR,Nst,wSTs_RHS,DATALOC) ; %

GEOproperties.RIGID_BODY_MODES = BasisUrb ;
GEOproperties.CENTROID = CENTROID ;
GEOproperties.INERTIA = INERTIA ;
GEOproperties.VOLUME = VOLUME ;


% --------------------BOUNDARY MATRICES
%disp('Computing shape function matrices for all boundary elements... (Nst)')
[ NelemB ,wSTb ] = ComputeNelemBoundALL(MESH.COOR,MESH.CNb,MESH.TypeElementB) ;
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
NstB_left = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngausB,MESH.CNb,nnode) ;
% Diagonal matrix with weights
ndimLOC = 1;
wDIAGb = CompWeightDiag(wSTb,ndimLOC)  ;
NstB_leftW = wDIAGb*NstB_left ;
GEOproperties.NstB_left = NstB_left;  % JAHO, 20-Nov-2023, see 
GEOproperties.wSTb = wSTb; 









% We are only interested in the contribution of such matrices at the faces
% defined in GID pre-process. Accordingly, we make
NstT_W_N_boundaries = cell(size(MESH.NODES_FACES)) ;

for iface = 1:length(NstT_W_N_boundaries)
    indlocal=  MESH.Indexes_faces_bnd_element{iface} ;
    indlocal_columns = small2large(indlocal,nnodeEb)  ;
    indlocal_rows =  small2large(indlocal,ngausB)  ;
    NstT_W_N_boundaries{iface} = NstB_leftW(indlocal_rows,:)'*NstB_right(indlocal_rows,indlocal_columns) ;
    
end




% Rigid body modes and area of each of the faces defined by the user in GID
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FACES = cell(size(MESH.Indexes_faces_bnd_element)) ;
for iface = 1:length(FACES)
    indELEMB = MESH.Indexes_faces_bnd_element{iface} ;
    CONNECTb_iface = MESH.CNb(indELEMB,:) ;
    [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(MESH.COOR,[],CONNECTb_iface,MESH.TypeElementB) ;
    nodes  = unique(CONNECTb_iface(:)) ;
    COOR_FACE = MESH.COOR(nodes,:) ; % Coordinates of this face
    COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
    %  BasisRrb = ConstructBasisRigidBody(COORrelA) ;
    FACES{iface}.AREA = AREA ;
    %   FACES{iface}.RIGID_BODY_MODES = BasisRrb ;
    FACES{iface}.COORrelA_global = COORrelA ;
    FACES{iface}.CENTROID = CentroidFA ;
    FACES{iface}.GeometricMassMatrix = Mst ; 
    
    [unitNORMALS,tTANGiniST] = NormalsAndTangentVectorsGaussPoints(MESH,CONNECTb_iface) ;
    
    FACES{iface}.UnitNormalAtGaussPoint = unitNORMALS ;
    FACES{iface}.UnitTangentAtGaussPoint = tTANGiniST ;
    
    
end
GEOproperties.FACES = FACES ;


