function [G,uBAR,DOFr,DOFm,AREA,R,DOFA,DOFB,M] = ...
    MIXED_CONDITION_faceCF(DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,DATA,NAME_COARSE)

if nargin == 0
    load('tmp2.mat')
end
Ub = [] ; 
% ---------
% FACE A   --- FINE MESH 
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

%%%%
%%%% FACE 2 - FINE MESH 
% ----------
iface=2 ;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COARSE MESH %
% ----------------------------------------------------------
[z,h,Mst_coarse,zNODES] = CoarseMeshVariable(NAME_COARSE,CentroidFA,COORrelA,nodesfA,COOR_FACE,AREA) ; 
% Slave DOFs 
DOFr = [DOFA(z); DOFA(h)] ;   
% Master DOFs 
DOFm = [DOFB(z); DOFB(h)] ; 


% --------------------------------------------
% Matrix Qbar
% -------------
% \Qbar \defeq \ident - \R_{z} \Par{ (\R_z^T \Mast \R_z)^{-1} \R_z^T \Mast}
% Mast = Mst_coarse
Mast = sparse(size(Mst_coarse,1)*ndim,size(Mst_coarse,2)*ndim) ;
for idim = 1:ndim
    Mast(idim:ndim:end,idim:ndim:end) = Mst_coarse ;
end
Ubar = Mast*R(z,:);
aB_coeff = (Ubar'*R(z,:))\Ubar';
Qbar = R(z,:)*aB_coeff ;
Qbar = speye(size(Qbar))-Qbar ;

% MATRIX H 
% ---------------------------
%  \H = \R_{h} \Par{ (\R_z^T \Mast \R_z)^{-1}   \R_z^T \Mast   }
H  = R(h,:)*aB_coeff ; 
% Matrix G 
% -------------------- 
% \G = \matcdos{\Qbar}{\zero}{-\H}{ \ident} 

G = sparse(length(DOFr),length(DOFm)) ;

iniROW = 1; 
finROW = iniROW + size(Qbar,1)-1 ; 
iniCOL = 1; 
finCOL = iniCOL + size(Qbar,2)-1 ; 
G(iniROW:finROW,iniCOL:finCOL) = Qbar ; 

% ---------- 
iniROW = finROW+1; 
finROW = iniROW + size(H,1)-1 ; 
iniCOL =  1; 
finCOL = iniCOL + size(H,2)-1 ; 
G(iniROW:finROW,iniCOL:finCOL) = -H ; 
% ----------
iniCOL =  finCOL+1; 
finCOL = length(DOFm); 
G(iniROW:finROW,iniCOL:finCOL) = speye(length(h),length(h)) ; 

uBAR  = zeros(length(DOFr),1) ; 





% 
% 
% 
%  QA = speye(length(DOFA),length(DOFA)) ;
% 
%  
%  Ubar = M*R;
%  coeff = (Ubar'*R)\Ubar';
%  Q = R*coeff ;
%  QB = speye(size(Q))-Q ;
%  
 
 

