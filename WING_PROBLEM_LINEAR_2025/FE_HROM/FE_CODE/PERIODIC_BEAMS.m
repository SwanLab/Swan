
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ; 
end
if iscell(a_B)
    a_B = cell2mat(a_B) ; 
end
da = a_A-a_B ;
%%%% FACE 1
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
CentoidFA = sum(COOR_FACE,1)/size(COOR_FACE,1); % Centroid
COORrelA = bsxfun(@minus,COOR_FACE',CentoidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ; % Basis matrix for rigid body motions, relative to centroid

% %%%% Identification of boundary nodes
% [CNloc1 ELEMloc] = ElemBnd(CONNECTb{1,iface},nodesfA) ;
% CN_BNDnodes = IdentifyBoundaryMesh(CNloc1) ;
% BNDnodes = unique(CN_BNDnodes) ; % Boundary nodes FACE 1
% % Global DOFs
% DOFA_b = small2large(BNDnodes,ndim) ;
% % Indices b (local)
% [dummy bindex] = intersect(DOFA,DOFA_b) ;
% [dummy iindex] = setdiff(1:length(DOFA),bindex) ;

%%%%%

% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired, 

%COOR_FACEB = COOR(nodesf,:) ;
%CentoidFB = sum(COOR_FACEB,1)/size(COOR_FACEB,1); % Center of gravity
%COORrelB = bsxfun(@minus,COOR_FACEB',CentoidFB')'; % Relative coordinates
%[IDX DISTANCES]= knnsearch(COORrelB,COORrelA) ;   % Matching nodes of face A with respect to B
%nodesfB = nodesf(IDX) ;  % Set of nodes of face B so that nodesfB(i) is paired with nodesfA(i)

DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B


%%%%%%55 NEW METHOD, Orthogonality edges of face 1
%--------------------------------------------------

% ORTHOGONALITY = 0 ;
% if ORTHOGONALITY ==0
M = [ R' , sparse(6,length(DOFA))
    speye(length(DOFA)) ,-speye(length(DOFA))
    ] ;

b = [ R'*(R*a_A)
    R*da
    ]  ;

% else
%     nconstraints = 6 +6 + length(DOFA) ;
%     MA = sparse(nconstraints,length(DOFA)) ;
%     MA(1:6,bindex) = R(bindex,:)' ;
%     MA(7:12,bindex) = R(bindex,:)' ;
%     MA(7:12,iindex) = R(iindex,:)' ;
%     MA(13:end,:) = speye(length(DOFA)) ;
%     
%     MB = sparse(nconstraints,length(DOFA)) ;
%     MB(13:end,:) = -speye(length(DOFA) );
%     
%     M = [MA, MB] ;
%     
%     b = zeros(size(MA,1),1) ;
%     b(1:6) = R(bindex,:)'*(R(bindex,:)*a_A) ;
%     b(7:12) = R'*(R*a_A) ;
%     b(13:end) = R*da ;
% end

%
%
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

G = -M_r\M_m ;
dR =  M_r\b ;


%%%%%%%%%%%%%%%%%%%%

%%% NEW METHOD ---------------




Gb = G ;









% %%%%%% ROTATION MATRIX of face B with respect centroid face A
% COOR_FACEB = COOR(nodesfB,:) ;
% COORrelAB = bsxfun(@minus,COOR_FACEB',CentoidFA')'; % Relative coordinates
% Rbar = ConstructBasisRigidBody(COORrelAB) ;

% %%%%%%%%%% Select 6 linearly independent columns from R_A
% [Xsub,c]=licols(R') ;
% %[Xsub,c]=licols(R') ;
% % Remaining DOFs (local numbering)
% f = 1:length(DOFA); f(c) = [] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Affine Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % IMPLEMENT 1
% % --------------------
% DOFr = [DOFA  ] ;  % Slave DOFs
% DOFm = DOFB(f) ;  % Master DOFs
% %
% % Matrix G
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 = -Rbar(c,:)'\Rbar(f,:)' ;
% %L3 =  -R(c,:)'\R(f,:)' ;
% G = sparse(length(DOFr),length(DOFm)) ;
% G(c,:) =  (L1) ;
% G(f,:) =  speye(length(f)) ;
%
% % Constant term
% m1 = -Rbar(c,:)'\(R'*(R*a_A)) ;
%
% dR = zeros(size(G,1),1) ;
% dR(c) = m1 ;
% dR = dR + R*da ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% imple = 1;
% if  imple==-1
%     %       DOFr = [DOFA ; DOFB(c)] ;  % Slave DOFs
%     DOFm = DOFB(f) ;  % Master DOFs
%     %
%     % Matrix G
%
%     GA = sparse(length(DOFA),length(DOFm)) ;
%     GA(f,:) = speye(length(f)) ;
%
%     G = [GA; sparse(length(c),length(DOFm))] ;
%     % Constant term
%     m1 = R(c,:)*a_A ;
%     m3 = R(c,:)*a_B ;
%
%     dA = zeros(length(DOFA),1) ;
%     dA(c) =  m1 ;
%     dA(f) = R(f,:)*da ;
%     dR3 = m3 ;
%     dR = [dA; dR3] ;
% elseif imple ==1
%     DOFr = [DOFA ; DOFB(c)] ;  % Slave DOFs
%     DOFm = DOFB(f) ;  % Master DOFs
%     %
%     % Matrix G
%
%     L1 = -Rbar(c,:)'\Rbar(f,:)' ;
%     L3 =  -R(c,:)'\R(f,:)' ;
%
%     GA = sparse(length(DOFA),length(DOFm)) ;
%     GA(c,:) = L1 ;
%     GA(f,:) = speye(length(f)) ;
%
%     G = [GA; sparse(L3)] ;
%     % Constant term
%     m1 = -Rbar(c,:)'\(R'*(R*a_A)) ;
%     m3 = R(c,:)'\(R'*(R*a_B)) ;
%
%     dA = zeros(length(DOFA),1) ;
% %     dA(c) = R(c,:)*da + m1 ;
% %     dA(f) = R(f,:)*da ;
% %     dR3 = m3 ;
% %     dR = [dA; dR3] ;
% elseif imple == 1

% else
%

%     DOFr = [DOFA ; DOFB(c)] ;  % Slave DOFs
%     DOFm = DOFB(f) ;  % Master DOFs
%     %
%     % Matrix G
%
%     L1 = -R(c,:)'\R(f,:)' ;
%     L3 =  -R(c,:)'\R(f,:)' ;
%
%     GA = sparse(length(DOFA),length(DOFm)) ;
%     GA(c,:) = L1 ;
%     GA(f,:) = speye(length(f)) ;
%
%     G = [GA; sparse(L3)] ;
%     % Constant term
%     m1 = R(c,:)'\(R'*(R*a_A)) ;
%     m3 = R(c,:)'\(R'*(R*a_B)) ;
%
%     dA = zeros(length(DOFA),1) ;
%     dA(c) = R(c,:)*da + m1 ;
%     dA(f) = R(f,:)*da ;
%     dR3 = m3 ;
%     dR = [dA; dR3] ;


%
%
% end





%%%%%%%%%%%%%%%%5
% M = [ R' , sparse(6,length(DOFA))
%       R' , sparse(6,length(DOFA))
%       speye(length(DOFA)) ,-speye(length(DOFA))
%       R'  , Rbar'] ;
%
% b = [ R'*(R*a_A)
%       R'*(R*a_B)
%       R*da
%        sparse(6,1)
%      ]  ;

% This option works (all linearly independent)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     M = [ R' , sparse(6,length(DOFA))
%         speye(length(DOFA)) ,-speye(length(DOFA))
%         ] ;
%
%     b = [ R'*(R*a_A)
%         R*da
%         ]  ;
%









%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%% ROTATION MATRIX of the structure
%     Centoid = sum(COOR,1)/size(COOR,1); % Centroid
%     COORrel = bsxfun(@minus,COOR',Centoid')'; % Coordinates relative to centroid
%     R_F = ConstructBasisRigidBody(COORrel) ; % Basis matrix for rigid body motions, relative to centroid
%
%
%     MODES_ROT_INCLUDE =[1,4,6] ;
%     M = [ R' , sparse(6,length(DOFA))
%           speye(length(DOFA)) ,-speye(length(DOFA))
%           R_F(DOFA,MODES_ROT_INCLUDE)'  , R_F(DOFB,MODES_ROT_INCLUDE)'] ;
%
%     b = [ R'*(R*a_A)
%           R*da
%          sparse(length(MODES_ROT_INCLUDE),1)]  ;
%






