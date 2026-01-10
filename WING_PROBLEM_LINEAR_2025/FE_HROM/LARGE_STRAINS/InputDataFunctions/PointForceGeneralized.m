function   Fpnt =  PointForceGeneralized(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ; 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/ImplemetationGeneralizedForces.mlx
% Generaliged forces on surfaces
% (patterned after /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/PointLoadsBeams.m)
% JAHO, 11-May-2022
if nargin == 0
  load('tmp1.mat')  
end

% PROPERTIES SURFACE UNDER CONSIDERATION 
nodesfA  = MESH.NODES_FACES{iface} ; 
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = MESH.COOR(nodesfA,:) ; % Coordinates of this face
CentroidFA = MESH.PROPERTIES_FACES{iface}.CENTROID  ; 
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;  % Rigid body modes 
Mst = MESH.PROPERTIES_FACES{iface}.GeometricMassMatrix  ; 
M = sparse(size(Mst,1)*ndim,size(Mst,2)*ndim) ;
for idim = 1:ndim
    M(idim:ndim:end,idim:ndim:end) = Mst ;
end

b_A_input = FORCE_INFO.AMPLITUDE(:) ; 

Fpnt = zeros(size(MESH.COOR,1)*ndim,1) ;
b_A = (R'*M*R)\b_A_input ;
% Nodal forces FACE A
Fpnt(DOFA) = M*R*b_A ;








 