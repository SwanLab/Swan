function [Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,...
    posgp,MESH,DATA,Kelas] = ...
    GeometricMatricesFun_EIFEbub(MESH,nstrain,PROPMAT,DATA)
% Matrix depending on the finite element mesh (i.e., that can be pre-computed)
% Adaptation of GeometricMatricesFun_EIFEbub.m  for the  case of
% CO-ROTATIONAL FORMULATION
% ----------------------------------------------------------------------------
% JAHO - 22-oct-2024, Sunday
if nargin == 0
    load('tmp1.mat')
end
%DATA = MESH.DATA;


if isfield(MESH,'COOR_SUPPORT')
    % Special type of finite element (uncoupled)
    % 25-April-2024
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
    CN_SUPPORT = MESH.CN_SUPPORT;
else
    CN_SUPPORT = [] ;
end

[Bmat_allelem,Nmat_allelem,WEIGHTSinteg,TRANSF_COORD,MESH.CN,DATA,Vrot,Kelem,MESH.CN_SUPPORT  ]=...
    B_N_matricesEIFEbub(MESH.COOR,MESH.CN, PROPMAT,MESH.MaterialType,DATA,MESH.TypeElement,CN_SUPPORT ) ;

% ELASTIC COARSE-SCALE STIFFNESS MATRIX
% 6-mARCH--2024
% sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
nelem = length(Kelem) ; % Number of coarse-scale elements
nnodeE_ext = size(DATA.MESHextended.CN,2); % Number of nodes per element, including the bubble node
ndim_ext= size(Kelem{1},1)/nnodeE_ext ;  % Number of DOFs per node (= number of bubble DOFs if this number is greater than ndim)
CN_ext = DATA.MESHextended.CN ; % extended connectivity matrix
nnode_ext = size(DATA.MESHextended.COOR,1) ;
disp(['Assemblying Kstiff global, using elemental Kelem'])
Kelas = AssemblyKGlobal(cell2mat(Kelem),nelem,nnodeE_ext,ndim_ext,CN_ext,nnode_ext) ;
Kelas = Kelas(DATA.MESHextended.DOFS_TO_KEEP,DATA.MESHextended.DOFS_TO_KEEP) ; % Eliminate Ghost DOFs
disp('...Done')
%   Kelas(:,DATA.MESHextended.DOFS_ghost)  =[] ; % Eliminate Ghost DOFs
%   Kelas(DATA.MESHextended.DOFS_ghost,:)  =[] ; % Eliminate Ghost DOFs



MESH.TRANSF_COORD = TRANSF_COORD;
MESH.Vrot = Vrot;
%1 )  Bmatrix. Matrix such that DeformationGradient
% ------------------------------------------------------------------------------------------------
if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1 && DATA.SMALL_STRAIN_KINEMATICS ==1
    wSTs = cell2mat(WEIGHTSinteg.INTforces) ;
    IDENTITY_F  = [] ;
    posgp = PROPMAT(1).EIFE_prop.INTforces.posgp' ;
    ngaus_STRESS = size(posgp,2) ;
    
    
    
    Bst_F = BstSmallStrains_EIFEbub(DATA.MESHextended,nstrain,Bmat_allelem) ;
    
else
    % JAHO, 25-MAY-2024
    wSTs = cell2mat(WEIGHTSinteg.INTforces) ;
    % IDENTITY_F  = [] ;
    posgp = PROPMAT(1).EIFE_prop.INTforces.posgp' ;
    ngaus_STRESS = size(posgp,2) ;
    
    [Bst_F,IDENTITY_F] = BstLargeStrains_EIFEbub(DATA.MESHextended,nstrain,Bmat_allelem,DATA) ;
end
MESH.DATA  = [] ;

% Stacked N-matrix. A matrix such that  u(x)  = Nst*d, where u(x) is the displacement at a given Gauss Point
wSTs_RHS = cell2mat(WEIGHTSinteg.BodyForces) ;
Nshapeelem  = cell2mat(Nmat_allelem) ;
posgp_RHS = PROPMAT(1).EIFE_prop.BodyForces.posgp' ;
ngaus_RHS = size(posgp_RHS,2) ;
Nst = AssemblyNGlobalBUB(Nshapeelem,DATA.MESHextended,ngaus_RHS,DATA ) ;

% COMPUTING THE CENTROID OF THE DOMAIN, AS WELL AS THE ROTATIONAL INERTIAS.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
DATALOC = [] ;
% DOFs BOUNDARY NODES

[BasisUrb,VOLUME,Rbar,CENTROID,INERTIA ]= ConstructBasisRigidBody_MASSM(MESH.COOR,Nst(:,DATA.MESHextended.DOFS_bLOC),wSTs_RHS,DATALOC) ; %

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
[nnode,ndim]  = size(MESH.COOR) ;
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


