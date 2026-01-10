%load('tmp1.mat')
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ; 
end
if iscell(a_B)
    a_B = cell2mat(a_B) ; 
end
%%%% FACE 1
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
% COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid
%%% Renumbering 
% Local coordinates  --> COOR_FACE 
% Local connectivities 
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
 
%%%%%

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

%%%%%%
indx = 1:3:length(DOFA) ;
indy = 2:3:length(DOFA) ;
indz = 3:3:length(DOFA) ;
  
%--------------------------------------------------
nconstraints = 6 +6 + 2*length(indx) ;
MA = sparse(nconstraints,length(DOFA)) ;
MA(1:6,indx) = R(indx,:)'*Mst ;
MA(1:6,indy) = R(indy,:)'*Mst;
MA(1:6,indz) = R(indz,:)'*Mst ;
MA(13:(13+length(indx)-1),indx) = speye(length(indx)) ;

MB = sparse(nconstraints,length(DOFA)) ;
MB(7:12,indx) = R(indx,:)' ;
MB(7:12,indy) = R(indy,:)';
MB(7:12,indz) = R(indz,:)';
MB((end-length(indx)+1):end,indx) = speye(length(indx)) ;
 
M = [MA, MB] ;

b = [Rorth*(R*a_A)
     Rorth*(R*a_B)
     R(indx,:)*a_A
      R(indx,:)*a_B] ; 


% [Xsub,IND]=licols(M') ;  % Not necessary !!!!!
% 
% % Linearly independent reactions --->IND
% M = M(IND,:);
% b = b(IND) ;

% Select invertible block
[Xsub,r]=licols(M) ; %
m = setdiff(1:size(M,2),r) ;
DOF_AB = [DOFA;DOFB] ;
DOFr = DOF_AB(r) ;
DOFm = DOF_AB(m) ;
M_r = M(:,r) ;
M_m = M(:,m) ;

Gb = -M_r\M_m ;
dR =  M_r\b ;

  