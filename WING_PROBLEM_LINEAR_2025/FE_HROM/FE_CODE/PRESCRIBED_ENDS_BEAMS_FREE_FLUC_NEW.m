
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1} ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2} ; % Rigid body amplitudes face 2
% Restrictions each end 
[r_A] = cellfun(@isempty,a_A) ;
[r_B] = cellfun(@isempty,a_B) ;
r_A = find(r_A == 0) ; 
r_B = find(r_B == 0) ; 
a_A = cell2mat(a_A(r_A)) ;
a_B = cell2mat(a_B(r_B)) ;
a_Ar = a_A(:) ; 
a_Br = a_B(:) ; 
% 
l_A = setdiff(1:6,r_A) ;  
l_B = setdiff(1:6,r_B) ; 





%%%% FACE 1
iface=1
nodesfA = unique(CONNECTb{iface}{1}) ;    % Nodes face1
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid 

%%%%%

% FACE 2
iface=2
nodesf = unique(CONNECTb{iface}{1}) ;    % Nodes face2
COOR_FACEB = COOR(nodesf,:) ;
CentoidFB = sum(COOR_FACEB,1)/size(COOR_FACEB,1); % Center of gravity
COORrelB = bsxfun(@minus,COOR_FACEB',CentoidFB')'; % Relative coordinates
[IDX DISTANCES]= knnsearch(COORrelB,COORrelA) ;   % Matching nodes of face A with respect to B
nodesfB = nodesf(IDX) ;  % Set of nodes of face B so that nodesfB(i) is paired with nodesfA(i)
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B


%%%%%%
indx = 1:3:length(DOFA) ;
indyz = setdiff(1:length(DOFA),indx ) ; 
%--------------------------------------------------
nconstraints = length(r_A) +length(r_B) + 2*length(indx) ;
MA = sparse(nconstraints,length(DOFA)) ;
MA(1:length(r_A),indx) = R(indx,r_A)' ;
MA(1:length(r_A),indyz) = R(indyz,r_A)';
sumR = length(r_A) + length(r_B) ; 
l = l_A ; 
if isempty(l)
    J_x = speye(length(indx))  ; 
    J_yz = sparse(length(indx),length(indyz));
else   
    H = R(indx,l)*((R(:,l)'*R(:,l))\R(:,l)') ;
    J_x = speye(length(indx))-H(:,indx) ; 
    J_yz =  -H(:,indyz) ; 
end
MA((sumR+1):(sumR+length(indx)),indx) = J_x;
MA((sumR+1):(sumR+length(indx)),indyz) = J_yz;


MB = sparse(nconstraints,length(DOFA)) ;
 

l = l_B ; 
if isempty(l)
    J_x = speye(length(indx))  ; 
    J_yz = sparse(length(indx),length(indyz));
else   
    H = R(indx,l)*((R(:,l)'*R(:,l))\R(:,l)') ;
    J_x = speye(length(indx))-H(:,indx) ; 
    J_yz =  -H(:,indyz) ; 
end
MB((end-length(indx)+1):end,indx) = J_x;
MB((end-length(indx)+1):end,indyz) = J_yz;
  

M = [MA, MB] ;

b = [R(:,r_A)'*(R(:,r_A)*a_Ar)
     R(:,r_B)'*(R(:,r_B)*a_Br)
     R(indx,r_A)*a_Ar
      R(indx,r_B)*a_Br] ; 


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

  