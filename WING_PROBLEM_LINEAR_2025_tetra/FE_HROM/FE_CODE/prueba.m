clc
clear all
load('ExampleBeamPeriodic_WS.mat')
load('tmp.mat')

a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ;
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ;
da = a_A-a_B ; 
%%%% FACE 1
iface=1
nodesfA = unique(CONNECTb{iface}{1}) ;    % Nodes face1
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid


% FACE 2
iface=2
nodesf = unique(CONNECTb{iface}{1}) ;    % Nodes face2
COOR_FACEB = COOR(nodesf,:) ;
CentoidFB = sum(COOR_FACEB,1)/size(COOR_FACEB,1); % Center of gravity
COORrelB = bsxfun(@minus,COOR_FACEB',CentoidFB')'; % Relative coordinates
[IDX DISTANCES]= knnsearch(COORrelB,COORrelA) ;   % Matching nodes of face A with respect to B
nodesfB = nodesf(IDX) ;  % Set of nodes of face B so that nodesfB(i) is paired with nodesfA(i)
DOFB = small2large(nodesf,ndim) ; % Set of DOFs face B

%%%%%% ROTATION MATRIX of face B with respect centroid face A
COOR_FACEB = COOR(nodesfB,:) ;
COORrelAB = bsxfun(@minus,COOR_FACEB',CentoidFA')'; % Relative coordinates
Rbar = ConstructBasisRigidBody(COORrelAB) ; 
 
%%%%%%%%%% Select 6 linearly independent columns from R_A
[Xsub,c]=licols(Rbar') ;
% Remaining DOFs (local numbering)
f = 1:length(DOFA); f(c) = [] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Affine Boundary conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DOFr = [DOFA ] ;  % Slave DOFs
DOFm = DOFB(f) ;  % Master DOFs 
%
% Matrix G 
L = -Rbar(c,:)'\Rbar(f,:)' ;
G = [sparse(L) ; speye(length(f))] ; 
% Constant term 
m = -Rbar(c,:)'\(R'*(R*a_A)) ; 
dR = zeros(size(G,1),1) ; 
dR(1:6) = m ; 
dR = dR + R*da ; 
   
% mAB =  (R_A(f,:)*a_A  -R_B(f,:)*a_B) ; 
% dR = [R_A(c,:)*a_A;     R_B(c,:)*a_B ;   mAB  ] ;
% % 
% % 
% Gb = [ sparse(6,length(f)) ;
%     sparse(6,length(f)) ;
%     speye(length(f))] ;


 