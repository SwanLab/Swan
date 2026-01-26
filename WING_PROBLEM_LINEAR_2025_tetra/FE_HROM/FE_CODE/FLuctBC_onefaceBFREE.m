function [DOFA,AREA,R,M,s,p,Ub,Ib] =  FLuctBC_onefaceBFREE(iface,DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,dBENDING,aA,DATA)

if nargin == 0
    load('tmp2.mat')
end

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

dBENDING = dBENDING(DOFA,:) ;
% 
% warning('Borarrrrarrrrr')
% load('tmp_3.mat','MODES')
% MODES = dBENDING ; 

% Basis Ub (bending fluctuating modes)
coeff = (R'*M*R)\(R'*M*dBENDING) ;
Ub = dBENDING -R*coeff ;
TOL = 1e-10 ;
if norm(Ub) < TOL
    error('No bending fluctuation modes. Perhaps Poisson s ratio is zero (or the tolerance is too low)')
end


PLOT_FLUCT_MODES = 0 ;

if PLOT_FLUCT_MODES == 1
    refFACE = 1;
    COORloc = COOR(DOMAINVAR.NODES_faces12{1,refFACE},:) ;
    CNref =  RenumberConnectivities( CONNECTb{refFACE},1:length(DOMAINVAR.NODES_faces12{1,refFACE}) );
    posgp = [] ;
    LEGENDG= ['FLUCT.'] ;
    NAME_MODES = [DATA.NameWS_bending_displacements(1:end-4),LEGENDG ];
    DATAMODES = [] ;
    GidPostProcessModes_dom(COORloc,CNref,TypeElementB,Ub,posgp,NAME_MODES,DATAMODES,LEGENDG);
end

ISSS = cellfun(@isempty,aA) ; 

if any(ISSS)
    
 error('There can be no free motion on face A')



end

UbBAR = M*Ub;
Rbar = M*R ; 
% Matrix Q
coeff = (UbBAR'*Ub)\UbBAR';
Q = Ub*coeff ;
Ib = speye(size(Q))-Q ;

%% Degrees of freedom p and s 
[~,p]=licols([Ub ]') ; %
s = setdiff(1:length(DOFA),p) ;





% % Subsets r and l
% [~,r]=licols(R(s,:)') ; %
% l = setdiff(1:length(s),r) ;
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