function [rDEF ]= AmplitudeReactions_PLATE(DATAROM,MESH2D,a,fextBEAMr,ndim,DATA_REFMESH_glo)

if nargin ==0
    load('tmp0z.mat')
end

%V = DATAROM.BasisInt ;  % Interface modes matrix
%ndim = size(V,2) ; % Number of entries for each node
nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
%nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)


nentities = length(DATAROM) ; % Number of joints/slices
KT = cell(1,nentities) ;
for ientities = 1:nentities
    % if iscell(DATAROM{ientities}.BasisInt)
    V = DATA_REFMESH_glo{ientities}.BasisINTfc   ;
    %else
    %   V1 = DATAROM{ientities}.BasisInt ;
    %  V2 = V1;
    %end
    %   T= cell(size(V)) ;
    %  for iface = 1:length(V)
    %T_1 = V1'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f1,:) ;
    T = V'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f,:) ;
    %KT{1,ientities} = DATAROM{ientities}.Kbeam*T_1' ;
    KT{ientities} = DATAROM{ientities}.Kbeam*T' ;
end

 

%
% VECTORIZED = 1;
%
%
% if  VECTORIZED == 0
%     rDEF = cell(nelem,1) ;  % Amplitude self-equilibrated reaction modes
%
%     for ielem = 1:nelem
%
%         ientity = MESH2D.MaterialType(ielem) ;
%
%         CNlocNOD = MESH2D.CN(ielem,:) ;
%         % CNloc = Nod2DOF(CNlocNOD,ndim) ;
%         rDEF{ielem} = -fextBEAMr{ielem}  ;
%         for inode = 1:nnodeE
%             NODE = CNlocNOD(inode) ;
%             DOFS = Nod2DOF(NODE,ndim) ;
%             rDEF{ielem} = rDEF{ielem} + KT{inode,ientity}*a(DOFS)  ;
%         end
%
%     end
%
% else
nmodes = size(fextBEAMr{1},1) ;   % number of modes reactions
fextBEAMr = cell2mat(fextBEAMr)  ;
rDEF = -fextBEAMr ;
rDEF  =reshape(rDEF,nmodes,[]) ;
for  ientity = 1:length(DATAROM)
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    CNlocNOD = MESH2D.CN(ELEMS,:) ;
    for inode = 1:size(CNlocNOD,2)
        NODE = CNlocNOD(:,inode) ;
        DOFS = Nod2DOF(NODE,ndim) ;
        aDOFS = a(DOFS) ;
        aDOFS = reshape(aDOFS,ndim,[]) ;
        COLUMNS_ini =  (inode-1)*ndim+1 ; 
         COLUMNS_fin = inode*ndim ;  
         COLUMNS = COLUMNS_ini:COLUMNS_fin ; 
        rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{ientity}(:,COLUMNS)*aDOFS  ;
    end
    
end




%
%
% end
%
%
% a;



