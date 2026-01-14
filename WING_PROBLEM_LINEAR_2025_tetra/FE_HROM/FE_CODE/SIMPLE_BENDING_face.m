function [DOFA,AREA,R,M,Ub,QA,QB] = ...
    SIMPLE_BENDING_face(DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,dBENDING,DATA,b_B_input)

if nargin == 0
    load('tmp.mat')
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
% FACE A 
% ---------
iface = 2; 
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 1

%dBENDING_old = dBENDING ;
dBENDING = [dBENDING(DOFA,:) ];
%dBENDING_B = [dBENDING_old(DOFB,:) ];

 iface = 1; 

if b_B_input(5) ~= 0   
    % Shear test in the y direction
    SELE = 1 ; 
    dBENDING  = dBENDING(:,1) ; 
    if any(b_B_input([2,6])) 
        error('Non-compatible values')
    end
    REMAINING_MODES= setdiff(1:6,[3,5]) ; 
    
elseif b_B_input(6) ~= 0   
     % Shear test in the z direction
     SELE = 2; 
    dBENDING  = dBENDING(:,2) ; 
    if any(b_B_input([3,5])) 
        error('Non-compatible values')
    end
    REMAINING_MODES= setdiff(1:6,[2,6]) ; 
end

% 
% warning('Borarrrrarrrrr')
% load('tmp_3.mat','MODES')
% MODES = dBENDING ; 

% Basis Ub (bending fluctuating modes)

coeff = (R'*M*R)\(R'*M*dBENDING) ;
Ub = dBENDING -R*coeff ;
TOL = 1e-10 ;
if norm(Ub) < TOL
    warning('No bending fluctuation modes. Perhaps Poisson s ratio is zero (or the tolerance is too low)')
Ub = [] ; 
end


PLOT_FLUCT_MODES = 1;

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


if isempty(Ub)
    QA = speye(length(DOFA),length(DOFA)) ;
else
UbBAR = M*Ub;
Rbar = M*R ; 
% ----------------------
% Matrix QA
% ----------------------
coeff = (UbBAR'*Ub)\UbBAR';
Q = Ub*coeff ;
QA = speye(size(Q))-Q ;

end

% ------ 
% Matrix QB 
INCLUDE_ALL_RIGID_BODY_MODES =1 ; 

if INCLUDE_ALL_RIGID_BODY_MODES ==1     
U = [Ub,R] ; 
else
    U = [Ub,R(:,REMAINING_MODES)] ; 
end


Ubar = M*U;
coeff = (Ubar'*U)\Ubar';
Q = U*coeff ;
QB = speye(size(Q))-Q ;
% 
% %% Degrees of freedom p and s 
% [~,p]=licols([Ub ]') ; %
% s = setdiff(1:length(DOFA),p) ;
% % Subsets r and l
% [~,r]=licols(R(s,:)') ; %
% l = setdiff(1:length(s),r) ;

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