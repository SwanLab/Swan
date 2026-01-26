function [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrix(COOR_FACE,nodesfA,CONNECTb,TypeElementB)

%COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% Geometric mass matrix
[Nst,wST] = GeometricMassMatrixSurface(CONNECTb,nodesfA,COOR_FACE,TypeElementB) ;  %
wSTdiag = CompWeightDiag(wST,1)  ;
Mst = (wSTdiag*Nst)'*Nst ;
% Recomputing centroid
CentroidFA = zeros(1,3) ;
AREA = sum(wST) ;
for idim = 1:3
    CentroidFA(idim) = wST'*(Nst*COOR_FACE(:,idim))/AREA ;
end