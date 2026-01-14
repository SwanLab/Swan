function [Nst,wSTb] = GeometricMassMatrixSurface(CONNECTb_face1,nodesf1,COOR_FACE,TypeElementB)
% Geometric mass matrix surface defined by connect.  CONNECTb_face1 and
% nodesf1, with coordinates COOR_FACE

% INPUTS 
% CONNECTb_face1  --> Connectivities of the face 
% nodesf1 --> List of nodes of the the face 
% COOR_FACE = COOR(nodesf1,:)


% Why to renumber connectivities ? 

[~, INPUT_NODES ]= sort(nodesf1) ;
CNface1= RenumberConnectivities(CONNECTb_face1  ,1:length(nodesf1)) ;
COOR_FACE = COOR_FACE(INPUT_NODES,:) ; 

[ Nelem wSTb ] = ComputeNelemBoundALL(COOR_FACE,CNface1,TypeElementB) ; 
% -------------------------------------------------------------------------
% Assembly 
nelem = size(CNface1,1) ; nnodeE = size(CNface1,2) ; ndim = 1; 
ngaus = size(Nelem,1)/nelem ; nnode = size(nodesf1,1) ;  
Nst = AssemblyNGlobal(Nelem,nelem,nnodeE,ndim,ngaus,CNface1,nnode) ;