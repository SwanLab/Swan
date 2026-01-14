function [G,uBAR,DOFr,DOFm,AREA,R] = PERIODIC_BEAMS_COARSE_FINE(DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    NameFileMeshCoarse,DATA)

if nargin == 0
    load('tmp1.mat')
end

% DATA = DefaultField(DATA,'FLUCTUATION_MODES_BENDING',[]) ;
% Ub = DATA.FLUCTUATION_MODES_BENDING ;

ndim = 3;
% da = a_A-a_B ;


%-------------------------------------------------------
% FINER MESH
%%%% FACE 1
iface=1 ;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
if  max(nodesfA-sort(nodesfA)) ~=0
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% Geometric mass matrix, centroid
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face

[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{1},TypeElementB) ;

%[CentroidFB,AREA,Mst] =CentroidGeometricMassMatrixNEW(DATA_REFMESH.COOR,nodesfB,DATA_REFMESH.CONNECTb{iface},TypeElementB) ;

%
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

% FACE 2
iface=2 ;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------
%   --------------
%   COARSE MESH
% ----------------

%New local varialble
MAKE_NODES_COINCIDENT = 1; % Moves the nodes of the coarse mesh so that the
% coincide with those of the fine mesh


SLICE.NAME = NameFileMeshCoarse ;
DATA3D = GeometrySlice(SLICE) ; %
% FACE 1
NODES_COARSE = DATA3D.NODES_FACES{1} ;
COOR_COARSE = DATA3D.COOR(NODES_COARSE,:) ; %

if  MAKE_NODES_COINCIDENT ==1
    % Finding set of nodes of COORrelA_coarse close to COORrelA
    COORrelA_coarse = bsxfun(@minus,COOR_COARSE',CentroidFA')'; % Coordinates relative to centroid
    [zNODES DIST]= knnsearch(COORrelA,COORrelA_coarse) ;
    
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
        
        error(['Refine the FINE mesh around the above list of nodes'])
    end
    
    
    COORrelA_coarse = COORrelA(zNODES,:) ;
    z = small2large(zNODES,ndim) ;
    COOR_COARSE = COOR_FACE(zNODES,:);
end
% Geometric mass matrix
[dummy, setBelemLOC]= ElemBnd(DATA3D.CNb,NODES_COARSE); % elements face "iface"
% Connectivities faces f1 and f2
CNb1= DATA3D.CNb(setBelemLOC,:) ;
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
    CentroidFA_coarse(idim) = wST_coarse'*(Nst_coarse*COOR_COARSE(:,idim))/AREA ;
end
COORrelA_coarse = bsxfun(@minus,COOR_COARSE',CentroidFA_coarse')'; % Coordinates relative to centroid

if abs(AREA_coarse-AREA) > 1e-2*AREA
    error('Fine and coarse meshes correspond to different geometries')
end

% Finding set of nodes of COORrelA_coarse close to COORrelA
if  MAKE_NODES_COINCIDENT ==0
    [zNODES DIST]= knnsearch(COORrelA,COORrelA_coarse) ;
    z = small2large(zNODES,ndim) ;
end

h = 1:size(COORrelA,1)*ndim;
h(z) = [] ;


% TAKE INTO ACCOUNT FLUCTUATION MODES 
% ------------------------------------
% 
% if ~isempty(Ub)
%     % Fluctuation modes (bending)
%     Ub = sparse(Ub) ; 
%     Uorth = sparse(zeros(size(Ub'))) ;
%     for idim =1:3
%         INDLOC =idim:3:size(Ub,1) ;
%         Uorth(:,INDLOC) = Ub(INDLOC,:)'*Mst ;
%     end
%     coeffs  = (Uorth*Ub)\Uorth ; 
%     Q = Ub*coeffs ; 
%     ident = speye(size(Q)) ; 
%     Rast = (ident-Q)\R ; 
%     
%     
% end



R_z = R(z,:) ;
%Rast_z = Rast(z,:) ; 
% Rbar --->
Rbar_z = zeros(size(R_z)) ;
for idim =1:ndim
    INDLOC =idim:ndim:size(R_z,1) ;
    Rbar_z(INDLOC,:) = Mst_coarse*R_z(INDLOC,:) ;
end


% Select 6 linearly independent rows from Rbar_z
[~,r]=licols(Rbar_z') ; %
l = setdiff(1:length(z),r) ;
%
J = inv(Rbar_z(r,:)')*Rbar_z(l,:)' ;
%b = inv(Rbar_z(r,:)')*Rbar_z'*R_z*a_A ;
DOFr = [DOFA(z(r)); DOFB(z(r)) ; DOFB(z(l)); DOFB(h)]  ;
DOFm = [DOFA(z(l)) ;DOFA(h) ];
G = sparse(length(DOFr),length(DOFm)) ;
uBAR = zeros(length(DOFr),1) ;

iini = 1;
ifin = size(J,1) ;
jini = 1 ;
jfin = size(J,2) ;
G(iini:ifin,jini:jfin) = - J ;
%uBAR(iini:ifin) = b  ;

iini = ifin+1;
ifin = iini+size(J,1)-1 ;
G(iini:ifin,jini:jfin) = + J ;
%uBAR(iini:ifin) = b - Rast_z(r,:)*da;

iini = ifin +1;
ifin = iini-1 + length(l) ;
G(iini:ifin,jini:jfin) = -speye(length(l));
%uBAR(iini:ifin) =   - Rast_z(l,:)*da;

iini = ifin +1;
jini = length(l) + 1;
G(iini:end,jini:end) = -speye(length(h));
%uBAR(iini:end) =   - Rast(h,:)*da;

