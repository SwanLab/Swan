function [G,uBAR,DOFr,DOFm,AREA,R] = PRESCRIBED_END_FREE_WARPING_fun(a_A,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0 
    load('tmp1.mat')
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


% [Nst,wST] = GeometricMassMatrixSurface(CONNECTb{1},nodesfA,COOR_FACE,TypeElementB) ; 
% wSTdiag = CompWeightDiag(wST,1)  ; 
% Mst = (wSTdiag*Nst)'*Nst ; 
% % Recomputing centroid 
% CentroidFA = zeros(1,3) ; 
% AREA = sum(wST) ; 
% for idim = 1:3 
%     CentroidFA(idim) = wST'*(Nst*COOR_FACE(:,idim))/AREA ; 
% end
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

 
% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired, 
 
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

% Rorth = R'*M --->  R(INDX)'*N'*W*N 

Rorth = zeros(size(R')) ; 
for idim =1:3 
    INDLOC =idim:3:size(R,1) ; 
    Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
end
Rbar = Rorth' ; 

 
 
nconstr = 6 ; 
% Select 6 linearly independent rows  
[~,r]=licols(Rorth) ; %
l = setdiff(1:length(DOFA),r) ;
% 
J = inv(Rbar(r,:)')*Rbar(l,:)' ; 
b = inv(Rbar(r,:)')*Rbar'*R*a_A ; 
DOFr = [DOFA(r)]  ; 
DOFm = DOFA(l) ; 
G = -J ; 
uBAR = b ; 
 