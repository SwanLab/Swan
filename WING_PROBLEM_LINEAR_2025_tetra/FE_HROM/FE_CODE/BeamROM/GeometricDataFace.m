function [R,GGG,CentroidFA,Mst,DOFA,Rwarping] = GeometricDataFace(iface,DATA_REFMESH,TypeElementB)

if nargin == 0
    load('tmp4.mat')
end

nodesfA = DATA_REFMESH.NODES_faces12{iface} ;  % Nodes face under consideration
COOR_FACE = DATA_REFMESH.COOR(nodesfA,:) ; % Coordinates of this face (global system)
ndim  =size(COOR_FACE,2 ) ; 

DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(DATA_REFMESH.COOR,nodesfA,DATA_REFMESH.CONNECTb{iface},TypeElementB) ; 
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid (global system)

% Now we have to rotate COORrelA to the system attatched to the face
if ~isempty(DATA_REFMESH.RotationMatrixFace{iface})
    COORrelA = (DATA_REFMESH.RotationMatrixFace{iface}'*COORrelA')' ;
end

R = ConstructBasisRigidBody(COORrelA) ;
% GEOMETRIC PROPERTIES 
nmodes = size(R,2) ; 
GGG = zeros(nmodes,nmodes)  ; 
for idim = 1:ndim
GGG = GGG + (R(idim:ndim:end,:)'*Mst*R(idim:ndim:end,:) ) ; 
end

% Let us compute the warping mode

 Rwarping = ConstructWarpingMode(COORrelA) ;
