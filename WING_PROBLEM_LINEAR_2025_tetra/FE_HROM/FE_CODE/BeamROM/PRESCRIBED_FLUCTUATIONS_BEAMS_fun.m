function [G,uBAR,DOFr,DOFm,AREA,R] = PRESCRIBED_FLUCTUATIONS_BEAMS_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0
    load('tmp.mat')
end

% DATA = DefaultField(DATA,'FLUCTUATION_MODES_BENDING',[]) ;
% Ub = DATA.FLUCTUATION_MODES_BENDING ;   % Not succesful...

ndim = 3;
%%%% FACE 1
iface=1 ;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
if  max(nodesfA-sort(nodesfA)) ~=0
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
% COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid
%%% Renumbering
% Local coordinates  --> COOR_FACE
% Local connectivities
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;


COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;


% FACE 2
iface=2 ;
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,

DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fluctuations
load(DATA.FLUCTUATIONS_BOUNDARY_NameWS,'BasisINTfluct')  ;
% Product of fluctuations times mass matrix
U  = BasisINTfluct ;
UbarT = zeros(size(BasisINTfluct')) ;
for idim =1:3
    INDLOC =idim:3:size(R,1) ;
    UbarT(:,INDLOC) = BasisINTfluct(INDLOC,:)'*Mst ;
end

% Check that UbarT*R = zero
CHECK = norm(UbarT*R) ;
TOL = 1e-10;
if CHECK > TOL
    error('Fluctuation modes are not orthogonal to rigid body modes')
end
% Matrix Q
coeff = (UbarT*U)\UbarT;
Q = U*coeff ;
Qbar = speye(size(Q))-Q ;



% Select 6 linearly independent rows from Ryz_n
[~,r]=licols(Qbar') ; %
l = setdiff(1:length(DOFA),r) ;

%
DOFr = [DOFA(r); DOFB(r) ]  ;
DOFm = [DOFA(l); DOFB(l)];
%
J = - Qbar(r,r)\Qbar(r,l) ;
uA =  Qbar(r,r)\(R(r,:)*a_A) ;
uB =  Qbar(r,r)\(R(r,:)*a_B) ;

G = sparse(length(DOFr),length(DOFm)) ;
uBAR = zeros(length(DOFr),1) ;

iini = 1;
ifin = size(J,1) ;
iiniC = 1;
ifinC = size(J,2) ;

G(iini:ifin,iiniC:ifinC) = J ;
uBAR(iini:ifin) = uA  ;

iini = ifin+1;
ifin = iini+size(J,1)-1 ;
iiniC = ifinC+1;
ifinC = iiniC+size(J,2)-1 ;

G(iini:ifin,iiniC:ifinC) = J ;
uBAR(iini:ifin) = uB ;


