function [Nst,wSTb] = GeometricMassMatrixSurfaceNEW(CONNECTb_face1,nodesf1,COOR_FACE,TypeElementB)
% Geometric mass matrix surface defined by connect.  CONNECTb_face1 and
% nodesf1, with coordinates COOR_FACE

% INPUTS 
% CONNECTb_face1  --> Connectivities of the face 
% nodesf1 --> List of nodes of the the face 

CNface1= RenumberConnectivities(CONNECTb_face1  ,1:length(nodesf1)) ;


[ Nelem wSTb ] = ComputeNelemBoundALL(COOR_FACE,CONNECTb_face1,TypeElementB) ; 
% -------------------------------------------------------------------------
% Assembly 
nelem = size(CONNECTb_face1,1) ; nnodeE = size(CONNECTb_face1,2) ; ndim = 1; 
ngaus = size(Nelem,1)/nelem ; nnode = size(nodesf1,1) ;  
Nst = AssemblyNGlobal(Nelem,nelem,nnodeE,ndim,ngaus,CNface1,nnode) ;