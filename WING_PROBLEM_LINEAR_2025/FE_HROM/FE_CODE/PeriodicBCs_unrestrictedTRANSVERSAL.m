function [G,dR,DOFr,DOFm] = PeriodicBCs_unrestrictedTRANSVERSAL(DOMAINVAR,COOR,CONNECTb,TypeElementB,a_A,a_B)

if nargin == 0
    load('tmp.mat')
end
%% Computing rotation matrix, mass matrix, and DOFS of face A and B
[R,DOFA,DOFB,Mst ] = DOFS_and_RotationMatrix_FACE_withMASS(DOMAINVAR,COOR,CONNECTb,TypeElementB) ; 
%%%%%%
indx = 1:3:length(DOFA) ;
indy = 2:3:length(DOFA) ;
indz = 3:3:length(DOFA) ; 
indyz = [indy,indz] ;
M = blkdiag(Mst,Mst) ; 
%--------------------------------------------------
c = [1,5,6] ; 
n = [2,3,4] ; 
nk = 4; 
R_xc = R(indx,c) ; 
R_yn = R(indy,n) ; 
R_zn = R(indz,n) ; 
R_ynk = R(indy,nk) ; 
R_znk = R(indz,nk) ;   
a_An = a_A(n) ; 
a_Bnk = a_B(nk) ; 
a_Ac = a_A(c) ; 
a_Bc = a_B(c) ; 
%----------------------------------------------------
nconstraints = 6 +1 + length(indx) ;
HA = sparse(nconstraints,length(DOFA)) ;
HA(1:3,indx) = R_xc'*Mst ;
HA(4:6,indy) = R_yn'*Mst;
HA(4:6,indz) = R_zn'*Mst;
ROWS = (nconstraints-length(indx) +1):nconstraints ; 
HA(ROWS,indx) = speye(length(indx)) ;

HB = sparse(nconstraints,length(DOFA)) ;
HB(7,indy) = R_ynk'*Mst ;
HB(7,indz) = R_znk'*Mst ;
ROWS = (nconstraints-length(indx) +1):nconstraints ; 
HB(ROWS,indx) =  -speye(length(indx)) ; 

H = [HA, HB] ;

b = zeros(nconstraints,1) ; 
b(1:3) = (R_xc'*Mst*R_xc)*a_Ac ; 
b(4:6) = (R_yn'*Mst*R_yn)*a_An + (R_zn'*Mst*R_zn)*a_An ; 
b(7) = (R_ynk'*Mst*R_ynk)*a_Bnk + (R_znk'*Mst*R_znk)*a_Bnk ; 
b(8:end) = R_xc*(a_Ac - a_Bc) ; 



 
%
%
% [Xsub,IND]=licols(M') ;  % Not necessary !!!!!
%
% % Linearly independent reactions --->IND
% M = M(IND,:);
%Z b = b(IND) ;

% Select invertible block
[Xsub,r]=licols(H) ; %
m = setdiff(1:size(H,2),r) ;
DOF_AB = [DOFA;DOFB] ;
DOFr = DOF_AB(r) ;
DOFm = DOF_AB(m) ;
H_r = H(:,r) ;
H_m = H(:,m) ;

G = -H_r\H_m ;
dR =  H_r\b ;


%%%%%%%%%%%%%%%%%%%%

%%% NEW METHOD ---------------


 








% %%%%%% ROTATION HATRIX of face B with respect centroid face A
% COOR_FACEB = COOR(nodesfB,:) ;
% COORrelAB = bsxfun(@minus,COOR_FACEB',CentoidFA')'; % Relative coordinates
% Rbar = ConstructBasisRigidBody(COORrelAB) ;

% %%%%%%%%%% Select 6 linearly independent columns from R_A
% [Xsub,c]=licols(R') ;
% %[Xsub,c]=licols(R') ;
% % ReHAining DOFs (local nuHBering)
% f = 1:length(DOFA); f(c) = [] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Affine Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % IMPLEMENT 1
% % --------------------
% DOFr = [DOFA  ] ;  % Slave DOFs
% DOFm = DOFB(f) ;  % HAster DOFs
% %
% % HAtrix G
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
%     DOFm = DOFB(f) ;  % HAster DOFs
%     %
%     % HAtrix G
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
%     DOFm = DOFB(f) ;  % HAster DOFs
%     %
%     % HAtrix G
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
%     DOFm = DOFB(f) ;  % HAster DOFs
%     %
%     % HAtrix G
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
%     %%%%% ROTATION HATRIX of the structure
%     Centoid = sum(COOR,1)/size(COOR,1); % Centroid
%     COORrel = bsxfun(@minus,COOR',Centoid')'; % Coordinates relative to centroid
%     R_F = ConstructBasisRigidBody(COORrel) ; % Basis HAtrix for rigid body motions, relative to centroid
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






