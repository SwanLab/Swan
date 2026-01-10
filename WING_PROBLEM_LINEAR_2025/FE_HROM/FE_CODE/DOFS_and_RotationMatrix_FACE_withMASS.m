function [R,DOFA,DOFB,Mst] = DOFS_and_RotationMatrix_FACE_withMASS(DOMAINVAR,COOR,CONNECTb,TypeElementB)

%%%% FACE 1
ndim = 3 ; 
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[Nst,wST] = GeometricMassMatrixSurface(CONNECTb{1},nodesfA,COOR_FACE,TypeElementB) ; 
wSTdiag = CompWeightDiag(wST,1)  ; 
Mst = (wSTdiag*Nst)'*Nst ; 
% Recomputing centroid 
CentroidFA = zeros(1,3) ; 
AREA = sum(wST) ; 
for idim = 1:3 
    CentroidFA(idim) = wST'*(Nst*COOR_FACE(:,idim))/AREA ; 
end
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;
 
% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,  
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B