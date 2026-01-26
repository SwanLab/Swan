function [rDEF ]= AmplitudeReactions_KINEM(DATAROM,MESH2D,a,fextBEAMr,ndim,DOFsKEEP)

if nargin ==0
    load('tmp2.mat')
end



%V = DATAROM.BasisInt ;  % Interface modes matrix
%ndim = size(V,2) ; % Number of entries for each node
nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)


nentities = length(DATAROM) ; % Number of joints/slices
KT = cell(1,nentities) ;





for ientities = 1:nentities    
    V = DATAROM{ientities}.BasisInt.BasisINTF ;
    [DOFsNOD,~]=  cellfun(@size,DATAROM{ientities}.BasisInt.BasisINTFall_cell) ;   
    T = V'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f,:) ;
    KT{ientities} = DATAROM{ientities}.Kbeam*T' ;    
end

nmodes = size(fextBEAMr{1},1) ;   % number of modes reactions
fextBEAMr = cell2mat(fextBEAMr)  ;
rDEF = -fextBEAMr ;
rDEF  =reshape(rDEF,nmodes,[]) ;

%%

if any(ndim - ndim(1))
    aALL = zeros(nnode*max(ndim),1) ;
    aALL(DOFsKEEP) = a;
    
    for  ientity = 1:length(DATAROM)
        ELEMS = find(MESH2D.MaterialType == ientity) ;
        CNlocNOD = MESH2D.CN(ELEMS,:) ;
        COLUMNS_ini  =1; 
        for inode = 1:nnodeE
            NODE = CNlocNOD(:,inode) ;
            DOFS = Nod2DOF(NODE,max(ndim)) ;
            aDOFS = aALL(DOFS) ;
            aDOFS = reshape(aDOFS,max(ndim),[]) ;
            nmodesLOC = size(KT{ientity},1) ;
            Ktloc = zeros(nmodesLOC,max(ndim)) ;
         %   COLUMNS_ini = iini + ndim(inode)-1 ; % (inode-1)*max(ndim)+1 ;
            COLUMNS_fin = COLUMNS_ini  + ndim(inode)-1  ; %   inode*max(ndim) ;
            COLUMNS = COLUMNS_ini:COLUMNS_fin ;
            Ktloc(1:nmodesLOC,1:ndim(inode)) = KT{ientity}(:,COLUMNS) ;
            COLUMNS_ini = COLUMNS_fin+1; 
            rDEF(:,ELEMS) =  rDEF(:,ELEMS) + Ktloc*aDOFS  ;
        end
    end
else
    ndim = ndim(1) ;
    for  ientity = 1:length(DATAROM)
        ELEMS = find(MESH2D.MaterialType == ientity) ;
        CNlocNOD = MESH2D.CN(ELEMS,:) ;
        for inode = 1:nnodeE
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
end
end



%
%
% %V = DATAROM.BasisInt ;  % Interface modes matrix
% %ndim = size(V,2) ; % Number of entries for each node
% nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
% nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
% %nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)
%
%
% nentities = length(DATAROM) ; % Number of joints/slices
% KT = cell(1,nentities) ;
% for ientities = 1:nentities
%     % if iscell(DATAROM{ientities}.BasisInt)
%     V = DATA_REFMESH_glo{ientities}.BasisINTfc   ;
%     %else
%     %   V1 = DATAROM{ientities}.BasisInt ;
%     %  V2 = V1;
%     %end
%     %   T= cell(size(V)) ;
%     %  for iface = 1:length(V)
%     %T_1 = V1'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f1,:) ;
%     T = V'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f,:) ;
%     %KT{1,ientities} = DATAROM{ientities}.Kbeam*T_1' ;
%     KT{ientities} = DATAROM{ientities}.Kbeam*T' ;
% end
%
%
%
% %
% % VECTORIZED = 1;
% %
% %
% % if  VECTORIZED == 0
% %     rDEF = cell(nelem,1) ;  % Amplitude self-equilibrated reaction modes
% %
% %     for ielem = 1:nelem
% %
% %         ientity = MESH2D.MaterialType(ielem) ;
% %
% %         CNlocNOD = MESH2D.CN(ielem,:) ;
% %         % CNloc = Nod2DOF(CNlocNOD,ndim) ;
% %         rDEF{ielem} = -fextBEAMr{ielem}  ;
% %         for inode = 1:nnodeE
% %             NODE = CNlocNOD(inode) ;
% %             DOFS = Nod2DOF(NODE,ndim) ;
% %             rDEF{ielem} = rDEF{ielem} + KT{inode,ientity}*a(DOFS)  ;
% %         end
% %
% %     end
% %
% % else
% nmodes = size(fextBEAMr{1},1) ;   % number of modes reactions
% fextBEAMr = cell2mat(fextBEAMr)  ;
% rDEF = -fextBEAMr ;
% rDEF  =reshape(rDEF,nmodes,[]) ;
% for  ientity = 1:length(DATAROM)
%     ELEMS = find(MESH2D.MaterialType == ientity) ;
%     CNlocNOD = MESH2D.CN(ELEMS,:) ;
%     for inode = 1:size(CNlocNOD,2)
%         NODE = CNlocNOD(:,inode) ;
%         DOFS = Nod2DOF(NODE,ndim) ;
%         aDOFS = a(DOFS) ;
%         aDOFS = reshape(aDOFS,ndim,[]) ;
%         COLUMNS_ini =  (inode-1)*ndim+1 ;
%          COLUMNS_fin = inode*ndim ;
%          COLUMNS = COLUMNS_ini:COLUMNS_fin ;
%         rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{ientity}(:,COLUMNS)*aDOFS  ;
%     end
%
% end
%
%
%
%
% %
% %
% % end
% %
% %
% % a;
%
%
%
