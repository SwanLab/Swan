function [DOFA,AREA,R,M,Ub,QA,QB,G,uBAR,DOFr,DOFm,DOFB] = ...
    MIXED_CONDITIONcoup_face(DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,INFOBENDING)

if nargin == 0
    load('tmp2.mat')
end
Ub = [] ;
% ---------
% FACE A
% ---------
iface = 1;
ndim = 3;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1


COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

M = sparse(size(Mst,1)*ndim,size(Mst,2)*ndim) ;
for idim = 1:ndim
    M(idim:ndim:end,idim:ndim:end) = Mst ;
end
QA = speye(length(DOFA),length(DOFA)) ;

%  % else
b  = INFOBENDING.b;

n = setdiff(1:6,b) ;
U = R(:,b) ;
%  % end
%

Ubar = M*U;
coeff = (Ubar'*U)\Ubar';
Q = U*coeff ;
QB = speye(size(Q))-Q ;

% ---------
% FACE B
% ---------
iface = 2;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rbar = M*R(:,n)  ; 
%% Degrees of freedom p and s
[~,p]=licols([Rbar ]') ; %
s = setdiff(1:length(DOFB),p) ;

% --------------------------------------------
DOFr = [DOFA; DOFB(p)]; 
DOFm = [DOFB(s)] ;
%%%%% 

J = Rbar(p,:)'\Rbar(s,:)' ; 
L  =QB(:,s) -QB(:,p)*J ; 
G = sparse([L; -J]) ; 
uBAR = zeros(length(DOFr),1) ; 


% % MATRIX A 
% A = [QA, -QB
%      R(:,n)'*M, zeros(2,size(QB,2))] ; %    zeros(4,size(QA,2)),  R(:,n)'*M] ; 
%   

% % Submatrices
% % -------------------------
% L_sp  = Ib(s,s)\Ib(s,p) ;
% % L_sr_p = L_sp(r,:) ;
% % L_sl_p = L_sp(l,:) ;
% hAB = Ib(s,s)\(R(s,:)*(aA-aB)) ;
% % hAB_r = hAB(r) ;
% % hAB_l = hAB(l) ;
%
% %
% J_rl = Rbar(s(r),:)'\Rbar(s(l),:)' ;
% J_rp = Rbar(s(r),:)'\Rbar(p,:)' ;
% b_As = Rbar(s(r),:)'\(Rbar'*R*aA) ;
% % MATR.s = s;
% % MATR.p = p ;
% % MATR.l = l ;
% % MATR.r = r;
%
%
%