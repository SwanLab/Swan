function [DOFr,DOFm,J,uA,AREA,R,Mst] =  FLuctuBC_oneface_FLZERO(iface,DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,aA)


ndim = 3;
 nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
 DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1

COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

%
DOFr = DOFA    ;  % All DOFs
DOFm = [];
%
J = [];
uA =  R*aA;
 
% % Matrix Q
% coeff = (UbarT*U)\UbarT;
% Q = U*coeff ;
% Qbar = speye(size(Q))-Q ;

% 
% 
% % Select n-6 linearly independent rows from Qbar
% SELECTION = 0;
% if  SELECTION == 1
%     [~,r]=licols([Qbar ]') ; %
%     l = setdiff(1:length(DOFA),r) ;
% else
%     [~,l]=licols([U ]') ; %
%     r = setdiff(1:length(DOFA),l) ;
% end

% %
% DOFr = [DOFA(r)  ]  ;
% DOFm = [DOFA(l) ];
% %
% J = - Qbar(r,r)\Qbar(r,l) ;
% uA =  Qbar(r,r)\(R(r,:)*aA) ;
