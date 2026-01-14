function [Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,...
    posgp,MESH] = ...
    GeometricMatricesFun_EIFE(MESH,nstrain,PROPMAT,DATA)
% Matrix depending on the finite element mesh (i.e., that can be pre-computed)
% Adaptation of GeometricMatricesFun.m  for the Empirical Interscale Finite
% Element Method 
% ----------------------------------------------------------------------------
% JAHO - 22-Nov-2020, Sunday
if nargin == 0
    load('tmp.mat')
end
%DATA = MESH.DATA;


[Bmat_allelem,Nmat_allelem,WEIGHTSinteg,TRANSF_COORD,MESH.CN,DATA,Vrot  ]=...
    B_N_matricesEIFE(MESH.COOR,MESH.CN, PROPMAT,MESH.MaterialType,DATA,MESH.TypeElement) ;


MESH.TRANSF_COORD = TRANSF_COORD; 
MESH.Vrot = Vrot; 
%1 )  Bmatrix. Matrix such that DeformationGradient
% ------------------------------------------------------------------------------------------------
if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1 && DATA.SMALL_STRAIN_KINEMATICS ==1
    wSTs = cell2mat(WEIGHTSinteg.INTforces) ; 
    IDENTITY_F  = [] ; 
    posgp = PROPMAT(1).EIFE_prop.INTforces.posgp' ; 
    ngaus_STRESS = size(posgp,2) ; 
    Bst_F = BstSmallStrains_EIFE(MESH,nstrain,Bmat_allelem) ;
else
    error('Option not implemented')
    [Bst_F,wSTs,ngaus_STRESS,IDENTITY_F,posgp] = BstLargeStrains(MESH,nstrain) ;
end
MESH.DATA  = [] ;

% Stacked N-matrix. A matrix such that  u(x)  = Nst*d, where u(x) is the displacement at a given Gauss Point
nnode = size(MESH.COOR,1); ndim = size(MESH.COOR,2); nelem = size(MESH.CN,1); nnodeE = size(MESH.CN,2) ;
% Number of Gauss points for integrating RHS and mass MATRIX may be different
% from number of Gauss points used in integrating stiffness matrix
% 
% MESH = DefaultField(MESH,'posgp_given',[]) ;
% if isempty(MESH.posgp_given)
%     TypeIntegrand = 'RHS';
% else
%     TypeIntegrand = {MESH.posgp_given,MESH.weights_given} ;
% end

%[~,~,shapef_RHS,~] = ComputeElementShapeFun(MESH.TypeElement,nnodeE,TypeIntegrand) ;

%[  wSTs_RHS,~, posgp_RHS] = ComputeW_RHS(MESH.COOR,MESH.CN,MESH.TypeElement,ndim,TypeIntegrand)  ;

wSTs_RHS = cell2mat(WEIGHTSinteg.BodyForces) ; 


%[ Nshapeelem,~  ] = ComputeNelemALL(MESH.TypeElement,nnodeE,ndim,nelem,TypeIntegrand) ;

Nshapeelem  = cell2mat(Nmat_allelem) ; 
posgp_RHS = PROPMAT(1).EIFE_prop.BodyForces.posgp' ; 

shapef_RHS = [] ; 
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


