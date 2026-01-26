function [Gb,dR,DOFr,DOFm,AREA,R] = SYMMETRIC_BEAMS_MASS_MATRIXfun(DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0 
    load('tmp1.mat')
end
 
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

% 
% if ~isempty(Ub)
%     % Fluctuation modes (bending)
%     Ub = sparse(Ub) ; 
%     Uorth = sparse(zeros(size(Ub'))) ;
%     %comprobar = zeros(2,6) ; 
%     for idim =1:3
%         INDLOC =idim:3:size(Ub,1) ;
%         Uorth(:,INDLOC) = Ub(INDLOC,:)'*Mst ;
%      %   comprobar = comprobar + Ub(INDLOC,:)'*Mst*R(INDLOC,:) ; 
%     end
%     coeffs  = (Uorth*Ub)\Uorth ; 
%     Q = Ub*coeffs ; 
%     ident = speye(size(Q)) ; 
%     Rast = (ident-Q)\R ; 
% else
%     Rast = R; 
%     
% end

Rast = R ; 

 
[Gb,dR,DOFr,DOFm] = BCs_BEAMS_SYMMETRIC(R,Rorth',Rast,DOFA,DOFB) ; 
 
% METHOD1 = 0 ;
% if METHOD1 ==1
%     % Select invertible block
%     [Xsub,r]=licols(M) ; %
%     %r = sort(r) ;
%     m = setdiff(1:size(M,2),r) ;
%     %m = sort(m) ;
%     DOF_AB = [DOFA;DOFB] ;
%     DOFr = DOF_AB(r) ;
%     DOFm = DOF_AB(m) ;
%     M_r = M(:,r) ;
%     M_m = M(:,m) ;
%     
%     G = -M_r\M_m ;
%     dR =  M_r\b ;
%     Gb = G ;
% 
%     
% else
%     
% end

 