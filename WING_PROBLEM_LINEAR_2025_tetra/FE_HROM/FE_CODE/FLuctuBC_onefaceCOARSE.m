function [DOFr,DOFm,G,uBAR,R] =  FLuctuBC_onefaceCOARSE(iface,DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,U,aA,R,Mst,NAMECOARSE,CENTROIDS)

if nargin == 0
    load('tmp.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINE  MESH
% -------------------------------------------------------
ndim = 3;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face, fine mesh

% Compute area FINE (for comparison purposes)
 [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;

%% Relative coordinates
% ---------------------
% Pseudo-centroid
CENTROID = sum(COOR_FACE',2)/size(COOR_FACE,1)' ; % (psesudo-)Centroid face1
% Which is therefore the actual centroid ?
% DIST = [] ;
% for icentr = 1:length(CENTROIDS)
%     DIST(icentr) = norm(CENTROIDS{icentr}(:)-pCENTROID(:)) ;
% end
% [~,IND] = min(DIST) ;
% CENTROID = CENTROIDS{IND} ;

COORrel = bsxfun(@minus,COOR_FACE',CENTROID)';


[z,h,Mst_coarse,zNODES,Mbar,AREA_coarse] = CoarseMeshVariableLOCAL(NAMECOARSE,COORrel,...
    nodesfA,iface,CENTROID,COOR_FACE) ;

%  if abs(AREA_coarse-AREA) > 1e-2*AREA
%      disp(['Fine and coarse meshes have different areas'])
%      disp(['This is a persistent error that should be subject to careful consideration'])
%      error('Fine and coarse meshes correspond to different geometries')
%  end

zNODES = nodesfA(zNODES) ;

% USE_IDENTITY_MATRIX = 1 ;
%
% if USE_IDENTITY_MATRIX == 1
% Mbar = speye(size(Mbar)) ;
% end
%%% MATRIX L
%  \L_A =  ( \U_{Az}^T \Mbar_A \U_{Az})^{-1}  \U_{Az}^T \Mbar_A
L  =(U(z,:)'*Mbar*U(z,:))\(U(z,:)'*Mbar) ;
% Matrix Qbar
%  \Qbar_A = \defeq \ident - \U_{Az} \L_A
Qbar = U(z,:)*L ;
Qbar = speye(size(Qbar))-Qbar ;
% Matrix J
% \J_A = -    \Qbar_{Arr}^{-1} \Qbar_{Arl}
% Select n-6 linearly independent rows from Qbar
[~,l]=licols(U(z,:)') ; %
r = setdiff(1:length(z),l) ;
J = -Qbar(r,r)\Qbar(r,l) ;
% Ind. term
% \u_{Azr} = \Qbar_{Arr}^{-1} \{\Qbar_{Ar}\R_{Az}} \a_A
uZr = Qbar(r,r)\(Qbar(r,:)*R(z,:)*aA) ;
% H
%\  \N_A \defeq (\L_{Ar} \J_A +  \L_{Al} )
N = L(:,r)*J + L(:,l) ;
%    \b \defeq \L_{Ar}\u_{Azr}  - \L_A \R_{Az} \a_A
b = L(:,r)*uZr - L*R(z,:)*aA ;
%  \H_A \=  \U_{Ah}  \N_A
H = U(h,:)*N ;
% \uBAR_{Ah} =  \R_{Ah} \a_A + \U_{Ah} \b
uH = R(h,:)*aA + U(h,:)*b ;

%%% SLAVE DOF
DOFr = [DOFA(z(r)); DOFA(h)] ;
% MASTER
DOFm = [DOFA(z(l))] ;

G = [J; H] ;
uBAR = [uZr;uH] ;


% %%% OTHER APPROACH
% OTHER_APPROACU =1 ;
% if  OTHER_APPROACU == 1
%     [~,l]=licols(U') ;
%     r = setdiff(1:length(DOFA),l) ;
% %
% %
% %     \item Choose just 6 points ($\l$), so that
% % \begin{equation}
% %  \d_{r}  = \R_{r} \a + \U_r \U_l^{-1} \d_l
% % \end{equation}
%
% DOFr = DOFA(r) ;
% DOFm = DOFA(l) ;
% G = U(r,:)*inv(U(l,:)) ;
% uBAR = R(r,:)*aA - U(r,:)*(inv(U(l,:))*(R(l,:)*aA)) ;
%
% end


end


function [z,h,Mst_coarse,zNODES,M,AREA_coarse] = CoarseMeshVariableLOCAL(NameFileMeshCoarse,COORrelA,...
    nodesfA,iface,CentroidFA,COOR_FACE)

%   --------------
%   COARSE MESH
% ----------------
ndim = 3;


SLICE.NAME = NameFileMeshCoarse ;
DATALOC.InterfacesDoNotMatch=1 ; 
DATA3D = GeometrySlice(SLICE,DATALOC) ; % % Mesh information, COARSE Mesh (COOR, CN,CNb)
% FACE iface
NODES_COARSE = DATA3D.NODES_FACES{iface} ;  % Nodes pertaining to "iface"
COOR_COARSE = DATA3D.COOR(NODES_COARSE,:) ; % Coordinates of the nodes
%[CentroidFA_coarse,AREA,~] =CentroidGeometricMassMatrixNEW(DATA3D.COOR,NODES_COARSE,DATA3D.CNb,DATA3D.TypeElementB) ;

% Finding set of nodes of COORrelA_coarse close to COORrelA
COORrelA_coarse = bsxfun(@minus,COOR_COARSE',CentroidFA)'; % Coordinates relative to centroid
[zNODES DIST]= knnsearch(COORrelA,COORrelA_coarse) ;
% zNODES --> List of nodes of the fine mesh (local numbering of the face)

[ZZZ ]= unique(zNODES) ;   % ZZZ = zNODES(iA) ; zNODES = ZZZ(iC)
if length(ZZZ) ~= length(zNODES)
    disp(['_________________________________________________________________________________________'])
    disp(['THERE IS ONE OR MORE NODES OF COARSE MESH PAIRED WITH THE SAME NODE OF THE FINE MESH'])
    disp(['_________________________________________________________________________________________'])
    
    % How to detect which node or nodes are repeated ?
    ZZZ = sort(zNODES) ;
    [Znr,iA,iC ]= unique(ZZZ) ;
    REPEATED = diff(iC) ;
    ZREPEAT = find(REPEATED==0) ;
    % This means that  ZREPEAT + 1 appears more than once in ZZZ
    RepeatedNodes = ZZZ(ZREPEAT) ;  % Numbering of face a
    GlobalNumbRepeated = nodesfA(RepeatedNodes) ;
    disp(['_________________________________________________________________________________________'])
    disp(['List of nodes: '])
    disp(['_________________________________________________________________________________________'])
    GlobalNumbRepeated
    clipboard('copy',num2str(GlobalNumbRepeated')) ;
    error(['Refine the FINE mesh around the above list of nodes'])
end

%----------------------------------------------------------------------
% Now we wish to construct the geometric mass matrix of the coarse mesh 
% ----------------------------------------------------------------------
%COORrelA_coarse = COORrelA(zNODES,:) ;  % Coordinates of the SELECTED NODES 
% first index of  zNODES 
z = small2large(zNODES,ndim) ;
COOR_COARSE = COOR_FACE(zNODES,:);  %  Coordinates of the SELECTED NODES 
% ---------------------------------------% ith row corresponds to point
% zNODES(i)
%%%%%%%%%%%%%%%%%% And which are the connectivities ? 
[dummy, setBelemLOC]= ElemBnd(DATA3D.CNb,NODES_COARSE); % elements face "iface"
% Connectivities face 
CNb1= DATA3D.CNb(setBelemLOC,:) ; 
% Geometric mass matrix

[Nst_coarse,wST_coarse] = GeometricMassMatrixSurface(CNb1,NODES_COARSE,COOR_COARSE,DATA3D.TypeElementB) ;  %
wSTdiag_coarse = CompWeightDiag(wST_coarse,1)  ;
Mst_coarse = (wSTdiag_coarse*Nst_coarse)'*Nst_coarse ;
CentroidFA_coarse = zeros(1,3) ;

NNNN = find((wST_coarse)<0) ;
if ~isempty(NNNN)
    disp('Unappropriate coarse mesh')
    
end

AREA_coarse = sum(wST_coarse) ;
for idim = 1:3
    CentroidFA_coarse(idim) = wST_coarse'*(Nst_coarse*COOR_COARSE(:,idim))/AREA_coarse ;
end
COORrelA_coarse = bsxfun(@minus,COOR_COARSE',CentroidFA_coarse')'; % Coordinates relative to centroid

% if abs(AREA_coarse-AREA) > 1e-2*AREA
%     error('Fine and coarse meshes correspond to different geometries')
% end

% % Finding set of nodes of COORrelA_coarse close to COORrelA
% if  MAKE_NODES_COINCIDENT ==0
%     [zNODES DIST]= knnsearch(COORrelA,COORrelA_coarse) ;
%     z = small2large(zNODES,ndim) ;
% end

h = 1:size(COORrelA,1)*ndim;
h(z) = [] ;

M = sparse(size(Mst_coarse,1)*ndim,size(Mst_coarse,1)*ndim) ;
for idim=1:ndim
    IND = idim:ndim:size(M,1) ;
    M(IND,IND) = Mst_coarse;
end


end
