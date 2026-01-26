function [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    PRESCRIBED_END_FREE_WARPING_funNEW(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0 
    load('tmp1.mat')
end


ndim = size(COOR,2); 

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
R = ConstructBasisRigidBody(COORrelA) ;

 
% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired, 
 
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

% Rorth = R'*M --->  R(INDX)'*N'*W*N 

Rorth = zeros(size(R')) ; 
for idim =1:ndim
    INDLOC =idim:ndim:size(R,1) ; 
    Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
end


error('Implement this function')
 
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
 