function [DOFr,DOFm,J,uA,AREA,R,Mst] =  FLuctuBC_oneface(iface,DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,BasisINTfluct,aA,BasisINTrb,M)
ndim = 3;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1

%USE_MASSMATRIX_SLICE = 1 ;
AREA = 0 ; 
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face


%if USE_MASSMATRIX_SLICE == 1
    Mst = M ;
    R = BasisINTrb ;
    
%else
    
%    [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;
%    COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
%    R = ConstructBasisRigidBody(COORrelA) ;
    
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fluctuations
% Product of fluctuations times mass matrix
U  = BasisINTfluct ;
UbarT = zeros(size(BasisINTfluct')) ;
for idim =1:ndim
    INDLOC =idim:3:size(R,1) ;
    UbarT(:,INDLOC) = BasisINTfluct(INDLOC,:)'*Mst ;
end

% Check that the SLICE Rigid body matrix coincide with the JOINT rigid body
% matrix
% if USE_MASSMATRIX_SLICE ==0
%     CHECK = norm(R-BasisINTrb)/norm(BasisINTrb) ;
%     CHECK_NORM = norm(UbarT*R);
%     TOL = 1e-8;
%     if CHECK > TOL
%         error('Fluctuation modes are not orthogonal to rigid body modes')
%     end
% end


% Matrix Q
coeff = (UbarT*U)\UbarT;
Qbar = U*coeff ;
Qbar = speye(size(Qbar))-Qbar ;



% Select n-6 linearly independent rows from Qbar
% SELECTION = 0;
% if  SELECTION == 1
%     [~,r]=licols([Qbar ]') ; %
%     l = setdiff(1:length(DOFA),r) ;
% else
    [~,l]=licols([U ]') ; %
    r = setdiff(1:length(DOFA),l) ;
%end

%
DOFr = [DOFA(r)  ]  ;
DOFm = [DOFA(l) ];
%
J = - Qbar(r,r)\Qbar(r,l) ;
uA =  Qbar(r,r)\(R(r,:)*aA) ;
