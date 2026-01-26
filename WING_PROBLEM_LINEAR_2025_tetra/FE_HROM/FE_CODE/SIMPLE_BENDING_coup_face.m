function [DOFA,AREA,R,M,Ub,QA,QB,G,uBAR,DOFr,DOFm,DOFB] = ...
    SIMPLE_BENDING_coup_face(DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,dBENDING,DATA,INFOBENDING,b_B_input)

if nargin == 0
    load('tmp2.mat')
end

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


% ---------
% FACE B
% ---------
iface = 2;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 1

%dBENDING_old = dBENDING ;
dBENDING = [INFOBENDING.dBENDING(DOFA,INFOBENDING.ORDER_MODE) ];
%dBENDING_B = [dBENDING_old(DOFB,:) ];
b= INFOBENDING.b ;
n = setdiff(1:6,b) ;



coeff = (R'*M*R)\(R'*M*dBENDING) ;
Ub = dBENDING -R*coeff ;
TOL = 1e-10 ;
if norm(Ub) < TOL
    warning('No bending fluctuation modes. Perhaps Poisson s ratio is zero (or the tolerance is too low)')
    Ub = [] ;
end


PLOT_FLUCT_MODES = 0;

if PLOT_FLUCT_MODES == 1
    %Ub =  dBENDING_B;
    refFACE = 1;
    COORloc = COOR(DOMAINVAR.NODES_faces12{1,refFACE},:) ;
    CNref =  RenumberConnectivities( CONNECTb{refFACE},1:length(DOMAINVAR.NODES_faces12{1,refFACE}) );
    posgp = [] ;
    LEGENDG= ['FLUCT.'] ;
    NAME_MODES = [DATA.NameWS_bending_displacements(1:end-4),LEGENDG ];
    DATAMODES = [] ;
    GidPostProcessModes_dom(COORloc,CNref,TypeElementB,Ub,posgp,NAME_MODES,DATAMODES,LEGENDG);
end


UbBAR = M*Ub;
Rbar = M*R ;
% ----------------------
% Matrix QA
% ----------------------
coeff = (UbBAR'*Ub)\UbBAR';
Q = Ub*coeff ;
QA = speye(size(Q))-Q ;

% ------
% Matrix QB s
U = [Ub,R(:,b)] ;
Ubar = M*U;
coeff = (Ubar'*U)\Ubar';
Q = U*coeff ;
QB = speye(size(Q))-Q ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h   = COOR(DOMAINVAR.NODES_faces12{end,2}(1),1) - COOR(DOMAINVAR.NODES_faces12{1,1}(1),1) ;
h = abs(h) ;
if b_B_input(5) ~= 0
    GAMMA = 1+ (b_B_input(3)/b_B_input(5))*h  ;
elseif b_B_input(6) ~= 0
    GAMMA = 1- (b_B_input(2)/b_B_input(6))*h  ;
else
    error('Option not implemented')
end


% MATRIX A
A = [QA,-QB
    R'*M, zeros(6,size(QB,2))
    zeros(length(n),size(QA,2)), R(:,n)'*M
    Ub'*M, -GAMMA*Ub'*M] ;



%%%% FACE 2
% ----------
iface=2 ;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 2
% ------------------------------------------------------------------


%% Degrees of freedom p and s
if isempty(Ub)
    p =  [] ;
elseif size(Ub,2) ==1
    [~, p] = max(abs(Ub)) ;
elseif size(Ub,2) >1
    [~,p]=licols([Ub ]') ;
end
s = setdiff(1:length(DOFA),p) ;


L = QA(s,s)\QA(s,p) ;
H = QA(s,s)\QB(s,:) ;
% MAtrix J
% --------------


% J = J_2/J_1
Ubar = M*Ub ;
J_1 = Ubar(p,:)'- Ubar(s,:)'*L;
J_2 = -Ubar(s,:)'*H + GAMMA*Ubar' ;
J = J_2/J_1 ;


% % Slave and master DOFs
DOFr = DOFA ;
DOFm = DOFB ;
G = sparse(length(DOFr),length(DOFm)) ;

G(s,:) = H-L*J ;
G(p,:) = J ;
uBAR =  zeros(length(DOFr),1) ;



