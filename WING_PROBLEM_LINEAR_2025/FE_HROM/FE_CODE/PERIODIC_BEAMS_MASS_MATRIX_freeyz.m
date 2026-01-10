
INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ; 
beta_p = INPUTS_LOC.FACTOR_PERIODICITY ; 
%load('tmp1.mat')
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
 
% FACE 2
iface=2 ; 
%nodesf = unique(CONNECTb{end,iface}) ;    % Nodes face2
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,  
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

% Rorth = R'*M --->  R(INDX)'*N'*W*N 

Rorth= zeros(size(R')) ; 
for idim =1:3 
    INDLOC =idim:3:size(R,1) ; 
    Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
end


%%%%%%
indx = 1:3:length(DOFA) ;
indy = 2:3:length(DOFA) ;
indz = 3:3:length(DOFA) ; 
%--------------------------------------------------
nconstraints = 6 +3 + length(indx) ;
MA = sparse(nconstraints,length(DOFA)) ;
MA(1:6,indx) = R(indx,:)'*Mst ;
MA(1:6,indy) = R(indy,:)'*Mst;
MA(1:6,indz) = R(indz,:)'*Mst;
MA(10:(10+length(indx)-1),indx) = beta_p*speye(length(indx)) ;

% 
MODESn = [2,3,4] ; 

MB = sparse(nconstraints,length(DOFA)) ;
MB(7:9,indx) = R(indx,MODESn)'*Mst ;
MB(7:9,indy) = R(indy,MODESn)'*Mst;
MB(7:9,indz) = R(indz,MODESn)'*Mst;
MB(10:(10+length(indx)-1),indx) = -speye(length(indx)) ;

M = [MA, MB] ;

b = [Rorth*(R*a_A)
    Rorth(MODESn,:)*(R*a_B)
    R(indx,:)*(a_A-a_B)] ;





% 
% M = [ Rorth , sparse(6,length(DOFA))
%     speye(length(DOFA)) ,-speye(length(DOFA))
%     ] ;
% 
% b = [ Rorth*(R*a_A)
%     R*da
%     ]  ;

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
%Z b = b(IND) ;

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






