function V = ModesCorners_Plates(DATA_REFMESH)

if nargin == 0
    load('tmp.mat')
end
% -----------------------------------------
icorner = 1; 
nodesCORNER= DATA_REFMESH.NODES_CORNERS{icorner} ;
COOR_CORNER = DATA_REFMESH.COOR(nodesCORNER,:) ; % Coordinates of this SIDE
zMID = 0.5*( max(COOR_CORNER(:,3))+ min(COOR_CORNER(:,3))) ; 
CentroidFA(1) = COOR_CORNER(1,1) ; 
CentroidFA(2) = COOR_CORNER(1,2) ; 
CentroidFA(3) = zMID; 
COORrelA = bsxfun(@minus,COOR_CORNER',CentroidFA')'; % Coordinates relative to centroid
% Rigid body modes (6)
R = ConstructBasisRigidBody(COORrelA) ;
V = R(:,1:5) ;  

end