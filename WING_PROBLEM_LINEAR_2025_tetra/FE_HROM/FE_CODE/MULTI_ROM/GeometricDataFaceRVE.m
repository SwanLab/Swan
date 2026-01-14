
function [R,GGG,CentroidFA,Mst] = GeometricDataFaceRVE(iface,DATA_REFMESH,TypeElementB)


ndim  =size(DATA_REFMESH.COOR,2) ;

nodesfA = DATA_REFMESH.NODES_faces{iface} ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = DATA_REFMESH.COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(DATA_REFMESH.COOR,nodesfA,DATA_REFMESH.CONNECTb{iface},TypeElementB) ;
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid

%%%%%%%%%%%
% MODF_change_anguloRVE. (22-MAY-2019)
% -----------------------------------
% Now we have to rotate COORrelA to the system attatched to the face
if ~isempty(DATA_REFMESH.RotationMatrixFace{iface})
    COORrelA = (DATA_REFMESH.RotationMatrixFace{iface}'*COORrelA')' ;
end

%%%%%%%%%%



R = ConstructBasisRigidBody(COORrelA) ;
%GEOMETRIC PROPERTIES
nrb = size(R,2) ;
GGG = zeros(nrb,nrb)  ;
for idim = 1:ndim
    GGG = GGG + (R(idim:ndim:end,:)'*Mst*R(idim:ndim:end,:) ) ;
end
 