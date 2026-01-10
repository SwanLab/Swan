function [rDEF ]= AmplitudeReactions_RVE(DATAROM,MESH2D,a,fextBEAMr,ndim,DOFsKEEP)

if nargin ==0
    load('tmp2.mat')
end


if isstruct(DATAROM{1}.BasisInt)
    % New method, Apr-2019 
    [rDEF ]= AmplitudeReactions_KINEM(DATAROM,MESH2D,a,fextBEAMr,ndim,DOFsKEEP) ; 
else
    
    
    
    
    %V = DATAROM.BasisInt ;  % Interface modes matrix
    %ndim = size(V,2) ; % Number of entries for each node
    nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
    nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
    nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)
    
    
    nentities = length(DATAROM) ; % Number of joints/slices
    KT = cell(nnodeE,nentities) ;
    
    
    
    
    
    for ientities = 1:nentities
        
        % if iscell(DATAROM{ientities}.BasisInt)
        V = DATAROM{ientities}.BasisInt ;
        %else
        %   V1 = DATAROM{ientities}.BasisInt ;
        %  V2 = V1;
        
        %end
     %  T= cell(size(V)) ;
        for iface = 1:length(V)
            %T_1 = V1'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f1,:) ;
           % T{iface} = V{iface}'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.fI{iface},:) ;
            %KT{1,ientities} = DATAROM{ientities}.Kbeam*T_1' ;
             KT{iface,ientities} = DATAROM{ientities}.Kbeam*(DATAROM{ientities}.Tcomp{iface})' ;
         %   KT{iface,ientities} = DATAROM{ientities}.Kbeam*T{iface}' ;
        end
        
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
    
    %%
    
    if any(ndim - ndim(1))
        aALL = zeros(nnode*max(ndim),1) ;
        aALL(DOFsKEEP) = a;
        
        for  ientity = 1:length(DATAROM)
            ELEMS = find(MESH2D.MaterialType == ientity) ;
            CNlocNOD = MESH2D.CN(ELEMS,:) ;
            for inode = 1:nnodeE
                NODE = CNlocNOD(:,inode) ;
                DOFS = Nod2DOF(NODE,max(ndim)) ;
                aDOFS = aALL(DOFS) ;
                aDOFS = reshape(aDOFS,max(ndim),[]) ;
                nmodesLOC = size(KT{inode,ientity},1) ;
                Ktloc = zeros(nmodesLOC,max(ndim)) ;
                Ktloc(1:nmodesLOC,1:ndim(inode)) = KT{inode,ientity} ;
                
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
                rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
            end
        end
    end
end

%
%
% end
%
%
% a;



