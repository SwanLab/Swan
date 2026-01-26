function [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    WARP_PERIODIC_BCSfun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

% COpy of PERIODIC_BEAMS_MASS_MATRIX for warping mode 15-Jan-2019
if nargin == 0 
    load('tmp1.mat')
end

a_A = [zeros(6,1); a_A]; 
a_B = [zeros(6,1); a_B]; 



ndim = size(COOR,2); 
da = a_A-a_B ;
%%%% FACE 1
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
if  max(nodesfA-sort(nodesfA)) ~=0 
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ; 
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid

BasisRrb = ConstructBasisRigidBody(COORrelA) ; % Rigid body modes 
Rwarping = ConstructWarpingMode(COORrelA) ;
%% Mass interface matrix 
M = sparse(size(Mst,1)*3,size(Mst,1)*3); 
for idim =1:ndim
    INDLOC =idim:ndim:size(Rwarping,1) ; 
    M(INDLOC,INDLOC) =  Mst ; 
end
% Make it orthogonal  
coeff = (BasisRrb'*M*BasisRrb)\(BasisRrb'*M*Rwarping) ; 
Rwarping = Rwarping - BasisRrb*coeff ; 
R = [BasisRrb,Rwarping] ; 
 
% FACE 2
iface=2 ; 
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,  
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% Rorth = R'*M --->  R(INDX)'*N'*W*N 
 

Rbar =M*R ; 

%%  \S = \Q \Mbar \R  = \Q \Rbar , where Q = diag(Q,Q,Q,Q)


% Matrix S (see )

[Gb,dR,DOFr,DOFm] = BCs_BEAMS_PERIODIC(DOFA,DOFB,R,a_A,a_B,Rbar) ; 
 
R = BasisRrb ; 
 