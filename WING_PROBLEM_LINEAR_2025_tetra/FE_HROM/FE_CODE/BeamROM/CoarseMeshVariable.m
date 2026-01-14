function [z,h,Mst_coarse,zNODES,M,AREA_coarse] = CoarseMeshVariable(NameFileMeshCoarse,CentroidFA,COORrelA,nodesfA,COOR_FACE,AREA)

%   --------------
%   COARSE MESH
% ----------------
ndim = 3;
%New local varialble



SLICE.NAME = NameFileMeshCoarse ;
DATA3D = GeometrySlice(SLICE) ; %
% FACE 1
NODES_COARSE = DATA3D.NODES_FACES{1} ;
COOR_COARSE = DATA3D.COOR(NODES_COARSE,:) ; %


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
    clipboard('copy',num2str(GlobalNumbRepeated')) ;
    error(['Refine the FINE mesh around the above list of nodes'])
end


COORrelA_coarse = COORrelA(zNODES,:) ;
z = small2large(zNODES,ndim) ;
COOR_COARSE = COOR_FACE(zNODES,:);

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
    CentroidFA_coarse(idim) = wST_coarse'*(Nst_coarse*COOR_COARSE(:,idim))/AREA_coarse ;
end
COORrelA_coarse = bsxfun(@minus,COOR_COARSE',CentroidFA_coarse')'; % Coordinates relative to centroid

if abs(AREA_coarse-AREA) > 1e-2*AREA
    error('Fine and coarse meshes correspond to different geometries')
end

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
