
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ;
end
if iscell(a_B)
    a_B = cell2mat(a_B) ;
end

%%%% FACE 1
iface=1
nodesfA = unique(CONNECTb{1,iface}) ;    % Nodes face1
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid

%%%%%

% FACE 2
iface=2
nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
COOR_FACEB = COOR(nodesf,:) ;
CentoidFB = sum(COOR_FACEB,1)/size(COOR_FACEB,1); % Center of gravity
COORrelB = bsxfun(@minus,COOR_FACEB',CentoidFB')'; % Relative coordinates
[IDX DISTANCES]= knnsearch(COORrelB,COORrelA) ;   % Matching nodes of face A with respect to B
nodesfB = nodesf(IDX) ;  % Set of nodes of face B so that nodesfB(i) is paired with nodesfA(i)
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B


%%%%%%


%--------------------------------------------------

EFFICIENT_IMPLEMENTATION = 1 ;

if EFFICIENT_IMPLEMENTATION == 1
    % GOAL --> Gb, dR 
    % 
   [Gb,dR,DOFr,DOFm] = BCs_BEAMS_freeFLUCyz(DOFA,DOFB,R,a_A,a_B)  ; 
 
else
    indx = 1:3:length(DOFA) ;

indyz = setdiff(1:length(DOFA),indx ) ;
    
    
    
    
    
    nconstraints = 6 +6 + 2*length(indx) ;
    MA = sparse(nconstraints,length(DOFA)) ;
    MA(1:6,indx) = R(indx,:)' ;
    MA(1:6,indyz) = R(indyz,:)';
    MA(13:(13+length(indx)-1),indx) = speye(length(indx)) ;
    
    MB = sparse(nconstraints,length(DOFA)) ;
    MB(7:12,indx) = R(indx,:)' ;
    MB(7:12,indyz) = R(indyz,:)';
    MB((end-length(indx)+1):end,indx) = speye(length(indx)) ;
    
    
    M = [MA, MB] ;
    
    b = [R'*(R*a_A)
        R'*(R*a_B)
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
    
    
end

